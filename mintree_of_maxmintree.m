function  mintree_out = mintree_of_maxmintree(tree, infield_name)
%  MINTREE_OF_MAXMINTREE computes the mintree of a maxtree or a mintree.
%
%  maxtree_out = MINTREE_OF_MAXMINTREE(tree, infield_name);
%
%  Input arguments:
%     tree:              	Structure for the maxtree or mintree
%     infield_name:         Field of the maxtree structure to analyze
%
%  Output argument:
%     maxtree_out:          Output mintree describing the connected
%                           components of infield_name in tree
%
%  EXAMPLE
%     maxtree_out = MINTREE_OF_MAXMINTREE(maxtree, 'Feature');
%
%  See also MAXTREE_OF_MAXMINTREE
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
if (~isfield(maxtree,infield_name))
    fprintf('Error in mintree_of_maxtree: Error in mintree_of_maxtree: The field "%s" does not exist.\n', infield_name);
    return;
end

% Process by duality
for node=1:length(tree)
    tree(node).(infield_name) = -1 * tree(node).(infield_name);
end

mintree_out = maxtree_of_maxmintree(tree, infield_name);

% Process by duality
for node=1:length(mintree_out)
    mintree_out(node).GrayLevel = -1 * mintree_out(node).GrayLevel;
end

end