function  tree_out = tree_Restitution(tree2, field2, tree1, field1)
%  TREE_RESTITUTION restitutes a tree field from a mintree or a maxtree 
%
%  tree_out = TREE_RESTITUTION(tree2, field2, tree1, field1);
%
%  Input arguments:
%     tree2:            Structure with the maxtree (or mintree) of the maxtree (or mintree) 
%     field2:           Name of the field of tree2 to restitute the tree1
%                       (likely to be 'GrayLevel')
%     tree1:            Structure with the maxtree (or mintree) defining 
%                       the connectivity for the restitution 
%     field1:           Name of the field of tree1 to store the restitution
%
%  Output argument:
%     tree_out:         Output tree1 with the added field1 
%
%  EXAMPLES
%       tree_out = TREE_RESTITUTION(tree2, 'GrayLevel', tree1, 'GrayLevel_Processed');
%
%  See also IMAGE_RESTITUTION
% 
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
tree_out = tree1;
if (~isfield(tree2,field2))
    fprintf('Error in tree_Restitution: The field "%s" does not exist in the tree passed as first argument.\n', field2);
    return;
end

%% Image restitution
for node2 = 1:length(tree2)
    for pix=1:tree2(node2).NumberOfPixels
        node1 = tree2(node2).Pixels(pix);
        tree_out(node1).(field1) = tree2(node2).(field2);
    end
end

end