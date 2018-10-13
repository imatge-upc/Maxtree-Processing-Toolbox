function  maxtree_out = maxtree_TDualReconstruct(maxtree, type, ref_field_name, mar_field_name, outfield_name)
%  MAXTREE_TDUALRECONSTRUCT computes the dual morphological Tree
%  reconstruction of mar_field_name with reference to ref_field_name in a
%  maxtree. The reconstruction process is either 'Root_to_Leaves' or
%  'Leaves_to_Root'. 
%
%  maxtree_out =  MAXTREE_TDUALRECONSTRUCT(maxtree, type, ref_field_name, ...
%                 mar_field_name, outfield_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     type:                 Type of reconstruction: 'Root_to_Leaves' or 
%                           'Leaves_to_Root'.
%     ref_field_name:       Name of the field to be used as reference 
%     mar_field_name:       Name of the field to be used as marker 
%     outfield_name:        Name of the output field  
%
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%       maxtree = MAXTREE_TDUALRECONSTRUCT(maxtree, 'Root_to_Leaves', ...
%                 'Graylevel', 'Marker', 'Out');
%
%  See also MAXTREE_TRECONSTRUCT
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
maxtree_out = maxtree;
if (~strcmp(type,'Root_to_Leaves')  && ...
    ~strcmp(type,'Leaves_to_Root')) 
    fprintf('Error in maxtree_TDualReconstruct: The reconstruction type "%s" does not exist.\n', type);
    return;
end
if (~isfield(maxtree, ref_field_name))
    fprintf('Error in maxtree_TDualReconstruct: The field "%s" does not exist.\n', ref_field_name);
    return;
end
if (~isfield(maxtree, mar_field_name))
    fprintf('Error in maxtree_TDualReconstruct: The field "%s" does not exist.\n', mar_field_name);
    return;
end

%% Compute the maximum of marker and reference values
ref = [maxtree.(ref_field_name)];
mar = [maxtree.(mar_field_name)];
tmp = num2cell(max(ref,mar));
[maxtree_out(:).(outfield_name)] = tmp{:};

if (strcmp(type,'Root_to_Leaves'))
    % Propagation from root to leaves
    for n=1:length(maxtree)
        % Propagate the value of node n to its children
        for i=1:length(maxtree_out(n).Children)
            node = maxtree_out(n).Children(i);
            val  = min(maxtree_out(n).(outfield_name),maxtree_out(node).(outfield_name));
            val  = max(val,maxtree_out(node).(ref_field_name));
            maxtree_out(node).(outfield_name) = val;
        end
    end
else
    % Propagation from leaves to root
    for n=length(maxtree):-1:2
        % Propagate the value of node n to its parent
        node = maxtree_out(n).Parent;
        val  = min(maxtree_out(n).(outfield_name),maxtree_out(node).(outfield_name));
        val  = max(val,maxtree_out(node).(ref_field_name));
        maxtree_out(node).(outfield_name) = val;
    end
end 

end