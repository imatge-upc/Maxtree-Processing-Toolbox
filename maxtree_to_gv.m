function  maxtree_to_gv(filename, layout, maxtree, varargin)
%  MAXTREE_TO_GV prints the Maxtree on a text file for vizualization with
%  graphviz. Graphviz is a package of open-source tools initiated by AT&T
%  Labs Research for drawing graphs specified in DOT language scripts. See:
%  http://www.graphviz.org
% 
%  MAXTREE_TO_GV(filename, layout, maxtree);
%
%  Input arguments:
%     filename:             Name fo the file to store the .gv file 
%     layout:               Graphviz layout (preferably 'dot', 'sfdp', 'neato' 
%                           or 'twopi')
%     maxtree:              Maxtree structure
%
%  Optional arguments:
%     'Subtree',val,        Only the subtree rooted in node "val" is displayed 
%                           (default 1)
%     'MaxNNodes',val,      Limit the number of nodes to be displayed to val
%                           (default all nodes are displayed)
%     'Attribute','name',   Name of the scalar field in the maxtree struct 
%                           to specify the gray level or color information 
%                           of the nodes (default 'GrayLevel')
%     'DisplayMode','name', Should be 'Color' or 'Gray' (default 'Gray') 
%     'Bounds',[min max],   Define the minimum and maximum values to be displayed 
% 
%  EXAMPLE
%       MAXTREE_TO_GV(filename, 'dot', maxtree);
%       MAXTREE_TO_GV(filename, 'sfdp', maxtree,'Subtree',15,'MaxNNodes',100,...
%           'Attribute','NumberOfPixels','DisplayMode','Color','Bounds',[0 10]);
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Process the variable list of input parameters
subtree = 1;
attribute_name  = 'GrayLevel';
display_mode    = 'Gray';
max_num_of_node = length(maxtree);
for i = 1 : 2 : length(varargin)
    switch varargin{i}
    case 'Subtree'
        subtree = varargin{i+1};
        if (subtree>length(maxtree)) 
            fprintf('Error in maxtree_to_gv: The subtree cannot be rooted in %d,', subtree);
            subtree = length(maxtree);
            fprintf(' using %d instead.\n', subtree);
        end
    case 'MaxNNodes'
        max_num_of_node = varargin{i+1};
    case 'Attribute'
        attribute_name = varargin{i+1};
        if (~isfield(maxtree,attribute_name)) 
            fprintf('Error in maxtree_to_gv: The Attribute "%s" does not exist \n', attribute_name);
            return;
        end
    case 'DisplayMode'
        display_mode = varargin{i+1};
        if (~(strcmp(display_mode,'Color')|| strcmp(display_mode,'Gray')) )
            fprintf('Error in maxtree_to_gv: The DisplayMode "%s" does not exist. It should be "Gray" or "Color".\n', display_mode);
            return;
        end
    case 'Bounds'
        bounds = varargin{i+1};
        min_attribute = bounds(1);
        max_attribute = bounds(2);
    end
end

%% Graphviz interface
fid = fopen(filename, 'w');
fprintf(fid, 'digraph G { \n');
fprintf(fid, 'graph [ordering= "out", layout="%s", rankdir="BT", nodesep=0.00  ] \n',layout);

%% Compute the max nad min values of the attribute
if (~exist('max_attribute','var'))
    max_attribute = maxtree(1).(attribute_name);
    min_attribute = maxtree(1).(attribute_name);
    for n=2:length(maxtree);
        max_attribute = max([max_attribute maxtree(n).(attribute_name)]);
        min_attribute = min([min_attribute maxtree(n).(attribute_name)]);
    end
end
norm_attribute = double(max_attribute)-double(min_attribute);

%% Compute the nodes to be displayed (depending on the subtree and the max number of nodes)
node_to_display = ones(length(maxtree),1);
if ((subtree~=1)||(max_num_of_node<length(maxtree)))    
    node_to_display(:)=0;
    fifo_queue = zeros(length(maxtree),1);
    fifo_start = 1; fifo_stop  = 1;
    fifo_queue(fifo_stop) = subtree; fifo_stop  = fifo_stop +1;
    i = 0;
    while((fifo_start < fifo_stop) && (i<max_num_of_node))
        node = fifo_queue(fifo_start); fifo_start = fifo_start +1;
        node_to_display(node) = 1; i = i+1; 
        for n=1:length(maxtree(node).Children)
            fifo_queue(fifo_stop) = maxtree(node).Children(n); fifo_stop  = fifo_stop +1;
        end
    end
end

%% Write the maxtree nodes
for n=1:length(maxtree);
    if (node_to_display(n)==1)
        if (strcmp(display_mode,'Gray'))
            str = dec2hex(uint8(250*(double(maxtree(n).(attribute_name)-double(min_attribute)))/norm_attribute),2);
            colorname = ['#' str str str];
            str = dec2hex(uint8(250*mod(0.5+(double(maxtree(n).(attribute_name)-double(min_attribute)))/norm_attribute,1)),2);
            colorname_font = ['#' str str str];
            %fprintf(fid, '%d [shape="circle", style="filled" , fillcolor="%s", color="%s", fontcolor="green" ];\n', maxtree(n).Node, colorname, colorname);
            fprintf(fid, '%d [shape="circle", style="filled" , fillcolor="%s", color="%s", fontcolor="%s" ];\n', maxtree(n).Node, colorname, colorname, colorname_font);
        else
            color = (1.0-(double(maxtree(n).(attribute_name))-double(min_attribute))/(double(max_attribute)-double(min_attribute)))*2/3;
            fprintf(fid, '%d [shape="circle", style="filled" , fillcolor="%f 1.00000 1.000000", color="%f 1.00000 1.000000"];\n', maxtree(n).Node, color, color);
        end
    end
end

%% write the edges
for n=1:length(maxtree);
    if (node_to_display(n)==1)
        for k=1:length(maxtree(n).Children)
            fprintf(fid, '%d -> %d [ color="#000000"];\n', maxtree(n).Node, maxtree(n).Children(k));
        end
    end
end

fprintf(fid, '}');
fclose(fid);
end