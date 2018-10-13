function  [branches, tree] = maxtree_to_branches(tree)
%  MAXTREE_TO_BRANCHES creates a "branches" map from a maxtree or mintree 
%  The branch map is an image where the pixel graylevel value indicates 
%  the node ID in the tree. The function will also change the order
%  of the children of each node so that they are in decreasing order of 
%  branch length. This allows a better visualization. The function output
%  is the branch map and the tree with the reordered children. 
% 
%  [branches, tree] = MAXTREE_TO_BRANCHES(tree);
%
%  Input arguments:
%     tree:              Structure with the maxtree or mintree 
%
%  Output arguments:
%     branches:          Branch map image
%     tree:              Input tree where children are ordered as a
%                        function of the branch length
%  EXAMPLE
%     [branches, tree] = MAXTREE_TO_BRANCHES(tree);
%
%  See also MAXTREE_BRANCHES_DISPLAY
%
%  Author: Sergi Liesegang Maria, Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Compute the distance to root for each node and the max branch length
tree(1).DistanceToRoot = 0;
number_of_branches     = 0;
for i = 2:length(tree)
    node  = tree(i).Node;
    nodeP = tree(i).Parent;
    tree(node).DistanceToRoot = tree(nodeP).DistanceToRoot + 1;
    if (isempty(tree(node).Children))
        number_of_branches = number_of_branches +1;
    end
end
height = 0;
for i = length(tree):-1:1
    node  = tree(i).Node;
    if (isempty(tree(node).Children))
        tree(node).MaxBranchLength = tree(node).DistanceToRoot;
    else
        tree(node).MaxBranchLength = 0;
        for n=1:length(tree(node).Children)
            nodeC = tree(node).Children(n);
            tree(node).MaxBranchLength = max(tree(node).MaxBranchLength, tree(nodeC).MaxBranchLength);
        end
    end
    height = max(height, tree(node).MaxBranchLength);
end
height = height + 1;

%% Re-oder the children as a function of their max branch length
for i = 1:length(tree)
    node  = tree(i).Node;
    if (~isempty(tree(node).Children))
        if (length(tree(node).Children)>1)
            nodeC = tree(node).Children(1);
            length_values = tree(nodeC).MaxBranchLength;
            for n=2:length(tree(node).Children)
                nodeC = tree(node).Children(n);
                length_values = [length_values tree(nodeC).MaxBranchLength];
            end
            [order_length_values, idx] = sort(length_values,'descend');
            tmp = tree(node).Children;
            tree(node).Children = [];
            for n=1:length(idx)
                tree(node).Children = [tree(node).Children tmp(idx(n))];
            end
        end
    end

end

%% Create the branch map
branches = int32(zeros(height,number_of_branches));
branch = 1;
lifo_queue    = zeros(length(tree),1); % queue with nodes id
lifo_queue(1) = 1;          % Put the root node in the queue
lifo_stop     = 1;
while(lifo_stop>0) 
    % Define the position of the node
    node = lifo_queue(lifo_stop); lifo_stop = lifo_stop - 1;
    branches(tree(node).DistanceToRoot+1,branch) = node;
    
    % Propagate the id from the left if we start a new branch
    if ((branch >1) && (branches(tree(node).DistanceToRoot,branch)==0))
        for i=1:tree(node).DistanceToRoot
            branches(i,branch) = branches(i,branch-1);
        end
    end
    
    % Put the children in the queue 
    for n=length(tree(node).Children):-1:1
        nodeC = tree(node).Children(n);
        lifo_stop = lifo_stop + 1;
        lifo_queue(lifo_stop) = nodeC;
    end
    % Increment the branch number if we process a leaf 
    if (isempty(tree(node).Children))
        branch = branch + 1;
    end
end

branches = flipud(branches);
end