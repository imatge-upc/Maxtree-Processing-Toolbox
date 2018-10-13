function  tree_out = maxtree_Prune(tree, attribute, threshold, decision)
%  MAXTREE_PRUNE prunes a mintree or a maxtree based on some attribute
%  values
%
%  tree_out = MAXTREE_PRUNE(tree, attribute, threshold, decision);
%
%  Input arguments:
%     maxtree:              Structure with the maxtree or mintree 
%     attribute:            Attribute on which the pruning is defined
%     threshold:            Value of the threshold on the attribute 
%     decision:             Type of decision: 'Direct, 'Min', 'Max' or
%                           'Subtract'
%  Output argument:
%     tree_out:             Output pruned tree 
%
%  EXAMPLE
%     maxtree_out = MAXTREE_PRUNE(maxtree, 'Area', 10, 'Direct');
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
tree_out = tree;
if (~strcmp(decision,'Min') && ...
    ~strcmp(decision,'Max') && ...
    ~strcmp(decision,'Direct') && ...
    ~strcmp(decision,'Subtract') ...
    ) 
    fprintf('Error in maxtree_Prune: The decision "%s" is not implemented.\n', decision);
    return;
end
if (~isfield(tree,attribute))
    fprintf('Error in maxtree_Prune: The field "%s" does not exist.\n', attribute);
    return;
end

%% Prune the tree "Direct"
if (strcmp(decision,'Direct') || strcmp(decision,'Subtract'))
     for node=length(tree):-1:1
         if tree_out(node).(attribute) <= threshold
            nodeP = tree_out(node).Parent;
            if (~isempty(nodeP))    % The root node cannot be removed
            
                % update Children
                for n=1:length(tree_out(node).Children)
                    nodeC = tree_out(node).Children(n);
                    tree_out(nodeC).Parent = nodeP;
                end
                tree_out(nodeP).Children(tree_out(nodeP).Children == node) = [];
                tree_out(nodeP).Children = [tree_out(nodeP).Children tree_out(node).Children];
                if (length(tree_out(nodeP).Children)==0)
                    tree_out(nodeP).Children=[];
                end

                % update NumberOfPixels
                tree_out(nodeP).NumberOfPixels = tree_out(nodeP).NumberOfPixels ...
                                                + tree_out(node).NumberOfPixels;
                % update Pixels
                tree_out(nodeP).Pixels = [tree_out(nodeP).Pixels tree_out(node).Pixels];

                tree_out(node).Node = [];
            end
         end     
     end
end

%% Modify the GrayLevel of tree_out in the case of the subtractive rule
if (strcmp(decision,'Subtract'))
    tree(1).Offset = 0;
    for node=2:length(tree)
        nodeP = tree(node).Parent;
        if isempty(tree_out(node).Node)
            tree(node).Offset = tree(nodeP).Offset ...
                - tree(node).GrayLevel + tree(nodeP).GrayLevel;
        else
            tree(node).Offset = tree(nodeP).Offset;
        end
    end
    for node=1:length(tree_out)
        tree_out(node).GrayLevel = tree_out(node).GrayLevel + tree(node).Offset;
    end
end

%% Prune the tree "Max"
if strcmp(decision,'Max')
     for node=length(tree):-1:1
         
         if ((tree_out(node).(attribute) <= threshold) && (length(tree_out(node).Children)==0))
            nodeP = tree_out(node).Parent;
            if (~isempty(nodeP))    % The root node cannot be removed
            
                % update Children
                for n=1:length(tree_out(node).Children)
                    nodeC = tree_out(node).Children(n);
                    tree_out(nodeC).Parent = nodeP;
                end
                tree_out(nodeP).Children(tree_out(nodeP).Children == node) = [];
                tree_out(nodeP).Children = [tree_out(nodeP).Children tree_out(node).Children];
                if (length(tree_out(nodeP).Children)==0)
                    tree_out(nodeP).Children=[];
                end

                % update NumberOfPixels
                tree_out(nodeP).NumberOfPixels = tree_out(nodeP).NumberOfPixels ...
                                                + tree_out(node).NumberOfPixels;
                % update Pixels
                tree_out(nodeP).Pixels = [tree_out(nodeP).Pixels tree_out(node).Pixels];

                tree_out(node).Node = [];
            end
         end     
     end
end

%% Prune the tree "Min"
if strcmp(decision,'Min')
     fifo_queue = zeros(length(tree),1); 
     for node=2:length(tree) % Start at 2 because the root node cannot be removed
         
         if ((tree_out(node).(attribute) <= threshold) && ...
                 (~isempty(tree_out(node).Node)))
            nodeP = tree_out(node).Parent;
            
            % Collapse the subtree hanging from node
            fifo_start = 1; fifo_stop  = 1;
            for n=1:length(tree(node).Children)
            	nodeC = tree_out(node).Children(n);
                fifo_queue(fifo_stop) = nodeC; 
                fifo_stop  = fifo_stop +1;
            end
            while(fifo_start < fifo_stop) 
                % Extract the node from the queue
                nodeL = fifo_queue(fifo_start); 
                fifo_start = fifo_start +1;
                
                % update NumberOfPixels
                tree_out(node).NumberOfPixels = tree_out(node).NumberOfPixels ...
                                            + tree_out(nodeL).NumberOfPixels;
                % update Pixels
                tree_out(node).Pixels = [tree_out(node).Pixels tree_out(nodeL).Pixels];

                tree_out(nodeL).Node = [];
                for n=1:length(tree(nodeL).Children)
                    nodeC = tree_out(nodeL).Children(n);
                    fifo_queue(fifo_stop) = nodeC; 
                    fifo_stop  = fifo_stop +1;
                end
            end
            
            % update Parent
            tree_out(nodeP).Children(tree_out(nodeP).Children == node) = [];
            if (length(tree_out(nodeP).Children)==0)
                tree_out(nodeP).Children=[];
            end
        
            % update NumberOfPixels
            tree_out(nodeP).NumberOfPixels = tree_out(nodeP).NumberOfPixels ...
                                            + tree_out(node).NumberOfPixels;
            % update Pixels
            tree_out(nodeP).Pixels = [tree_out(nodeP).Pixels tree_out(node).Pixels];

            tree_out(node).Node = [];
         end     
     end
end


%% Clean (remove void nodes)
current_node = 1;
LUT  = zeros(length(tree_out),1);
LUTi = zeros(length(tree_out),1);
for node=1:length(tree_out)
    if (~isempty(tree_out(node).Node))
        LUT (current_node) = node;
        LUTi(node) = current_node;
        current_node = current_node +1;
    end
end
for node=1:current_node-1;
    tree_tmp(node)      = tree_out(LUT(node));
    tree_tmp(node).Node = node;
    if (~isempty(tree_out(LUT(node)).Parent))
        tree_tmp(node).Parent = LUTi(tree_out(LUT(node)).Parent);
    end
    tree_tmp(node).Children = [];
    for nc=1:length(tree_out(LUT(node)).Children)
        tree_tmp(node).Children = [tree_tmp(node).Children LUTi(tree_out(LUT(node)).Children(nc))];
    end
end
       
tree_out = tree_tmp;

end