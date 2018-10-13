function  maxtree_out = maxtree_RegExtrema(maxtree, extrema, infield_name, outfield_name)
%  MAXTREE_REGEXTREMA identifies regional extrema of attributes on maxtrees
%  or mintrees
%
%  maxtree_out = MAXTREE_REGEXTREMA(maxtree, extrema, infield_name, outfield_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     extrema:              'Min' or 'Max'
%     infield_name:         Field of the maxtree structure on which to
%                           compute the extremas
%     outfield_name:        Field to define the position of the extremas
%                           Min/Max = 1, otherwise 0
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE
%     maxtree_out = MAXTREE_REGEXTREMA(maxtree, 'Min', 'GrayLevel', 'Marker');
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
maxtree_out = maxtree;
if (~strcmp(extrema,'Min') && ~strcmp(extrema,'Max')) 
    fprintf('Error in maxtree_RegExtrema: The regiona extrema cannot be "%s". It should be "Min" or "Max".\n', extrema);
    return;
end
if (~isfield(maxtree,infield_name))
    fprintf('Error in maxtree_RegExtrema: The field "%s" does not exist.\n', infield_name);
    return;
end

% Minima are found by duality
if (strcmp(extrema,'Min'))
    tmp = num2cell(-1*[maxtree.(infield_name)]);
    [maxtree(:).(infield_name)] = tmp{:};
end

%% Propagation from root to leaves
maxtree_out(1).(outfield_name) = 1;
for n=1:length(maxtree)
    % Compare the node n value to the values of its children
    for i=1:length(maxtree(n).Children)
        node = maxtree(n).Children(i);
        if (maxtree(n).(infield_name) < maxtree(node).(infield_name))
            maxtree_out(n).(outfield_name)    = 0;
            maxtree_out(node).(outfield_name) = 1;
        elseif (maxtree(n).(infield_name) > maxtree(node).(infield_name))
            maxtree_out(node).(outfield_name) = 0;
        end
    end
    % The label of the node n may have changed, so propagate the final value
    for i=1:length(maxtree(n).Children)
        node = maxtree(n).Children(i);
        if (maxtree(n).(infield_name) == maxtree(node).(infield_name))
            maxtree_out(node).(outfield_name) = maxtree_out(n).(outfield_name);
        end
    end
end

%% Propagation from leaves to root
for n=length(maxtree):-1:2
    % Compare the value of the node n and the value of its parent
    node = maxtree(n).Parent;
    if ((maxtree_out(node).(outfield_name)==1) && ...
        (maxtree_out(n).(outfield_name)==0) && ...
         (maxtree(n).(infield_name) >= maxtree(node).(infield_name)))
        maxtree_out(node).(outfield_name) = 0;
    end
end

%% Second propagation from root to leaves 
% (Necessary as the propagation from leaves to root may have changed the values of nodes)
for n=1:length(maxtree)
    for i=1:length(maxtree(n).Children)
        node = maxtree(n).Children(i);
        if (maxtree(n).(infield_name) == maxtree(node).(infield_name))
            maxtree_out(node).(outfield_name) = maxtree_out(n).(outfield_name);
        end
    end
end


end