function  maxtree_out = maxtree_of_maxmintree(tree, infield_name)
%  MAXTREE_OF_MAXMINTREE computes the maxtree of a maxtree or a mintree
%
%  maxtree_out = MAXTREE_OF_MAXMINTREE(tree, infield_name);
%
%  Input arguments:
%     tree:              	Structure for the maxtree or mintree
%     infield_name:         Field of the maxtree structure to analyze
%
%  Output argument:
%     maxtree_out:          Output maxtree describing the connected
%                           components of infield_name in tree
%
%  EXAMPLE
%     maxtree_out = MAXTREE_OF_MAXMINTREE(maxtree, 'Feature');
%
%  See also MINTREE_OF_MAXMINTREE
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
[tree(:).Labeled]=deal(0);
if (~isfield(tree,infield_name))
    fprintf('Error in maxtree_of_maxmintree: The field "%s" does not exist.\n', infield_name);
    return;
end

%% Sort the values of the initial tree
values = [tree(:).(infield_name)];
[sorted_values, sorted_index] = sort(values,'descend');
unique_values = fliplr(unique(sorted_values));

%% Flooding
fifo_queue  = zeros(length(tree),1);
fifo_start  = 1; fifo_stop  = 1;
value_start = 1;
value_end   = 1;
while (sorted_values(value_end)==sorted_values(value_start))
    value_end = value_end +1;
    if (value_end > length(sorted_values)) 
        break
    end
end
current_node = 0;

for ii = 1:length(unique_values)
    % Proceed with the flooding of the highest gray levels
    for n=value_start:value_end-1
        node = sorted_index(n);
        if (tree(node).Labeled==0) % We have found a node without label
        
            % Check whether the node is connected to a branch via its parent
            % or children 
            new_node = 0;
            nodeP = tree(node).Parent;
            if (~isempty(nodeP))
                if (tree(nodeP).Labeled~=0)
                    fifo_queue(fifo_stop) = node; 
                    fifo_stop  = fifo_stop +1;
                    new_node = 1;
                end
            end
            if (new_node == 0)
                for nc=1:length(tree(node).Children)
                    nodeC = tree(node).Children(nc);
                    if (tree(nodeC).Labeled~=0)
                        fifo_queue(fifo_stop) = node; 
                        fifo_stop  = fifo_stop +1;
                        new_node = 1;
                        break;
                    end
                end
            end

            % The unlabeled node is already connected to a branch, so
            % process it
            if (new_node ==1)             
                % A new node has to be created 
                current_node = current_node +1;
                maxtree_out(current_node).Node = current_node;
                maxtree_out(current_node).GrayLevel = sorted_values(n);
                maxtree_out(current_node).Parent = [];
                maxtree_out(current_node).Children = [];
                maxtree_out(current_node).Pixels = [];

                % Process its flatzone
                neighbor_list = [];
                while(fifo_start < fifo_stop)
                    node = fifo_queue(fifo_start);
                    fifo_start = fifo_start +1;
                    tree(node).Labeled = current_node;
                    maxtree_out(current_node).Pixels = [maxtree_out(current_node).Pixels node];

                    % Propagation to Parent nodes
                    nodeP = tree(node).Parent;
                    if (~isempty(nodeP))
                        if ((tree(nodeP).(infield_name)==tree(node).(infield_name))&& ...
                            (tree(nodeP).Labeled==0))
                            fifo_queue(fifo_stop) = nodeP; 
                            fifo_stop  = fifo_stop +1;
                        end
                        if (tree(nodeP).Labeled~=0)
                            % Look for the first (closest) ancestor
                            nodeL = tree(nodeP).Labeled;
                            while(~isempty(maxtree_out(nodeL).Parent))
                                nodeL = maxtree_out(nodeL).Parent;
                            end
                            if (maxtree_out(nodeL).GrayLevel~=maxtree_out(current_node).GrayLevel)
                                maxtree_out(current_node).Children = [maxtree_out(current_node).Children nodeL];
                                maxtree_out(nodeL).Parent = current_node;
                            end
                            neighbor_list = [neighbor_list nodeL];
                        end
                    end

                    % Propagation to Children nodes
                    for nc=1:length(tree(node).Children)
                        nodeC = tree(node).Children(nc);
                        if ((tree(nodeC).(infield_name)==tree(node).(infield_name))&& ...
                                (tree(nodeC).Labeled==0))
                            fifo_queue(fifo_stop) = nodeC; 
                            fifo_stop  = fifo_stop +1;
                        end
                        if (tree(nodeC).Labeled~=0)
                            % Look for the first (closest) ancestor
                            nodeL = tree(nodeC).Labeled;
                            while(~isempty(maxtree_out(nodeL).Parent))
                                nodeL = maxtree_out(nodeL).Parent;
                            end
                            if (maxtree_out(nodeL).GrayLevel~=maxtree_out(current_node).GrayLevel)
                                maxtree_out(current_node).Children = [maxtree_out(current_node).Children nodeL];
                                maxtree_out(nodeL).Parent = current_node;
                            end
                            neighbor_list = [neighbor_list nodeL];
                        end

                    end
                end
                % do we have to merge the newly created CC with an already
                % existing one?  
                neighbor_list=unique(neighbor_list);
                for n1=1:length(neighbor_list)
                    nodeL = neighbor_list(n1);
                    if ((nodeL~=current_node))
                        if(maxtree_out(current_node).GrayLevel==maxtree_out(nodeL).GrayLevel)
                            Label = tree(maxtree_out(nodeL).Pixels(1)).Labeled;
                            maxtree_out(nodeL).Children = [maxtree_out(nodeL).Children maxtree_out(current_node).Children];
                            maxtree_out(nodeL).Pixels   = [maxtree_out(nodeL).Pixels   maxtree_out(current_node).Pixels];
                            for m=1:length(maxtree_out(current_node).Children)
                                node1 = maxtree_out(current_node).Children(m);
                                maxtree_out(node1).Parent = nodeL;
                            end
                            for m=1:length(maxtree_out(current_node).Pixels)
                                node1 = maxtree_out(current_node).Pixels(m);
                                tree(node1).Labeled = Label;
                            end
                            maxtree_out(current_node)   = [];
                        end
                    end
                end
            end   
        end
    end
    

    % Process the new connected components at the current gray level 
    for n=value_start:value_end-1
        node = sorted_index(n);
        if (tree(node).Labeled==0)
            fifo_queue(fifo_stop) = node; 
            fifo_stop  = fifo_stop +1;
            
            % A new node has to be created 
            current_node = current_node +1;
            maxtree_out(current_node).Node = current_node;
            maxtree_out(current_node).GrayLevel = sorted_values(n);
            maxtree_out(current_node).Parent = [];
            maxtree_out(current_node).Children = [];
            maxtree_out(current_node).Pixels = [];
            
            % Process its flatzone
            while(fifo_start < fifo_stop)
                node = fifo_queue(fifo_start);
                fifo_start = fifo_start +1;
                tree(node).Labeled = current_node;
                maxtree_out(current_node).Pixels = [maxtree_out(current_node).Pixels node];
                
                % Propagation to Parent nodes
                nodeP = tree(node).Parent;
                if (~isempty(nodeP))
                    if ((tree(nodeP).(infield_name)==tree(node).(infield_name))&& ...
                        (tree(nodeP).Labeled==0))
                        fifo_queue(fifo_stop) = nodeP; 
                        fifo_stop  = fifo_stop +1;
                    end
                end
                
                % Propagation to Children nodes
                for nc=1:length(tree(node).Children)
                    nodeC = tree(node).Children(nc);
                    if ((tree(nodeC).(infield_name)==tree(node).(infield_name))&& ...
                            (tree(nodeC).Labeled==0))
                        fifo_queue(fifo_stop) = nodeC; 
                        fifo_stop  = fifo_stop +1;
                    end
                end
            end
        end
    end

    % Define the range of the next gray level value to process
    value_start = value_end;
    if (value_start<=length(sorted_values))
        while (sorted_values(value_end)==sorted_values(value_start))
            value_end = value_end +1;
            if(value_end>length(sorted_values))
                break;
            end
        end
    end
end

%% Debug
% for node=1:length(tree)
%     node1 = tree(node).Labeled;
%     list1 = maxtree_out(node1).Pixels;
%     if (max(list1==node)==0)
%         error=1;
%     end
% end

%% Clean (remove void nodes) and reorder the tree (so that the root is the first node)
current_node = 1;
LUT  = zeros(length(maxtree_out),1);
LUTi = zeros(length(maxtree_out),1);
for node=length(maxtree_out):-1:1
    if (~isempty(maxtree_out(node).Node))
        LUT (current_node) = node;
        LUTi(node) = current_node;
        current_node = current_node +1;
    end
end
for node=1:current_node-1;
    maxtree_tmp(node)      = maxtree_out(LUT(node));
    maxtree_tmp(node).Node = node;
    if (~isempty(maxtree_out(LUT(node)).Parent))
        maxtree_tmp(node).Parent = LUTi(maxtree_out(LUT(node)).Parent);
    end
    maxtree_tmp(node).Children = [];
    for nc=1:length(maxtree_out(LUT(node)).Children)
        maxtree_tmp(node).Children = [maxtree_tmp(node).Children LUTi(maxtree_out(LUT(node)).Children(nc))];
    end
end
% Define the NumberOfPixels field
for node=1:length(maxtree_tmp)
    maxtree_tmp(node).NumberOfPixels = length(maxtree_tmp(node).Pixels);
end
maxtree_out = maxtree_tmp;

end