function  maxtree_plot(maxtree, node, field_name, varargin)
%  MAXTREE_PLOT plots the field_name values of a Maxtree from the node
%  "node" to the root.
%
%  MAXTREE_PLOT(maxtree, node, field_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     node:                 Id of the node to start the plot (the plot goes
%                           up to the root node)
%     field_name:           Field to be displayed
%     varargin:             Plot option of the type '-r', '--b', etc. 
%
%  EXAMPLE
%     MAXTREE_PLOT(maxtree, 10, 'GrayLevel', '-r');
%
%  See also MAXTREE_TO_BRANCHES, MAXTREE_BRANCHES_DISPLAY
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
if (node > length(maxtree)) 
    fprintf('Error in maxtree_plot: The node "%d" does not exist in the tree.\n', node);
    return;
end
if (~isfield(maxtree, field_name))
    fprintf('Error in maxtree_plot: The field "%s" does not exist.\n', field_name);
    return;
end


%% Propagation from the node to root
val = maxtree(node).(field_name);
x = cellstr(sprintf('%d',node));
while (node>1)    
    node = maxtree(node).Parent;
    val = [val maxtree(node).(field_name)];
    x = [x; cellstr(sprintf('%d',node))];
end


if (isempty(varargin))
    plot(val);
else
    plot(val,varargin{1});
end

title(field_name);
step = int32(max(1,length(x)/20));
set(gca,'XTick',1:step:length(x))
set(gca,'XTickLabel',x(1:step:length(x)))
xlabel('Node ID');
ylabel('Value');

end