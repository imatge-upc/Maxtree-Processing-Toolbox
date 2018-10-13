function  image = image_Restitution(tree, field_name, width, height)
%  IMAGE_RESTITUTION restores an image from a mintree or a maxtree.
%
%  image = IMAGE_RESTITUTION(tree, field, width, height);
%
%  Input arguments:
%     tree:              Structure with the maxtree or mintree 
%     field_name:        name of the field to create the image
%     width:             width of the image
%     height:            height of the image
%
%  Output argument:
%     image:             Output image 
% 
%  EXAMPLE 
%     ima = IMAGE_RESTITUTION(tree, 'GrayLevel', 256, 256);
%
%  See also TREE_RESTITUTION
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
if (~isfield(tree,field_name))
    fprintf('Error in image_Restitution: The field "%s" does not exist.\n', field_name);
    return;
end

%% Image restitution
image = zeros(width, height);
for node = 1:length(tree)
    for pix=1:tree(node).NumberOfPixels
        image(tree(node).Pixels(pix)+1) = tree(node).(field_name);
    end
end

end