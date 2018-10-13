function  image = maxtree_branches_display(maxtree, branches, field_name, varargin)
%  MAXTREE_BRANCHES_DISPLAY displays a maxtree field as an image following
%  the branch mapping.
%
%  image = MAXTREE_BRANCHES_DISPLAY(maxtree, branches, field_name);
%
%  Input arguments:
%     maxtree:              Maxtree or mintree structure
%     branches:             Branches mapping created by maxtree_to_branches
%     field_name:           Field to be displayed
%  Optional arguments:
%     'AboveTreeValue',val, Value of the upper part of the image that does
%                           not correspond to actual nodes (default 0)
%  Output argument:        
%     image:                Image describing the evolution of field_name
%                           values along the tree branches
%
%  EXAMPLE 
%     image = MAXTREE_BRANCHES_DISPLAY(maxtree, branches, field_name);
%
%  See also MAXTREE_TO_BRANCHES
%
%  Authors: Sergi Liesegang Maria, Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
if (~isfield(maxtree, field_name))
    fprintf('Error in maxtree_branches_display: The field "%s" does not exist.\n', field_name);
    return;
end
AboveTreeValue = 0;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
    case 'AboveTreeValue'
        AboveTreeValue = varargin{i+1};
    end
end


%% Image creation
image      = zeros(size(branches));
image(:,:) = AboveTreeValue;
for i=1:size(image,1)
    for j=1:size(image,2)
        if(branches(i,j)>0)
            image(i,j) = maxtree(branches(i,j)).(field_name);
        end
    end
end


end