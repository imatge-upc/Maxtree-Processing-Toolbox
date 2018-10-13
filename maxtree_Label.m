function  maxtree_out = maxtree_Label(maxtree, infield_name, outfield_name)
%  MAXTREE_LABEL assigns a different label to all connected components that
%  are different from 0.
%
%  maxtree_out = MAXTREE_LABEL(maxtree, infield_name, outfield_name);
%
%   Input arguments:
%     maxtree:              Maxtree structure
%     infield_name:         Field of the maxtree structure to label
%     outfield_name:        Field to output the labels
%
%   Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%       maxtree = MAXTREE_LABEL(maxtree, infield_name, outfield_name));
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
maxtree_out = maxtree;
if (~isfield(maxtree,infield_name))
    fprintf('Error in maxtree_Label: The field "%s" does not exist.\n', infield_name);
    return;
end

%% Define the label of the root
tmp = num2cell(int16(zeros(1,length(maxtree))));
[maxtree_out.(outfield_name)]=tmp{:};

if ((maxtree(1).(infield_name))==0)
    Label_num = 0;
else
    Label_num = 1;
end
maxtree_out(1).(outfield_name)=Label_num;

%% Propagation from root to leaves
for n=2:length(maxtree)
    if ((maxtree(n).(infield_name)) ~= 0)
        node = maxtree(n).Parent;
        if ((maxtree(node).(infield_name)) == 0)
            Label_num = Label_num + 1;
            maxtree_out(n).(outfield_name)=Label_num;
        else
            maxtree_out(n).(outfield_name)=maxtree_out(node).(outfield_name);
        end
    end
    
end

end