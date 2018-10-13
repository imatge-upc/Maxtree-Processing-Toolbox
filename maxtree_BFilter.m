function  maxtree_out = maxtree_BFilter(maxtree, filter, parameter, aggregation, infield_name, outfield_name)
%  MAXTREE_BFILTER computes a Branch Filter on attributes of Maxtrees. In
%  Branch filters, each branch is individually filtered (estimation step)
%  and then the resulting estimated values are aggregated.
%
%  maxtree_out = MAXTREE_BFILTER(maxtree, filter, parameter, aggregation,...
%                infield_name, outfield_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     filter:               Type of (estimation) filters, possible values are: 
%                           'Erosion','Dilation', 'Mean','Median'
%     parameter:            Size of the neighborhood in the branch
%     aggregation:          Type of aggregation, possible values are:
%                           'Mean','Median', 'Max', 'Min''
%     infield_name:         Name of the field defining the input data
%     outfield_name:        Name of the field defining the output data
%
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%     maxtree = MAXTREE_BFILTER(maxtree,'Opening',2,'Median','GrayLevel','Out');
%
%  See also MAXTREE_GFILTER, MAXTREE_TFILTER, MAXTREE_BFILTER2
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
maxtree_out = maxtree;
if (~strcmp(filter,'Erosion')  && ...
    ~strcmp(filter,'Dilation') && ...  
    ~strcmp(filter,'Mean')  && ...
    ~strcmp(filter,'Median')...
    ) 
    fprintf('Error in maxtree_BFilter: The filter "%s" does not exist.\n', filter);
    return;
end
if (~strcmp(aggregation,'Min')  && ...
    ~strcmp(aggregation,'Max') && ...  
    ~strcmp(aggregation,'Mean')  && ...
    ~strcmp(aggregation,'Median')...
    ) 
    fprintf('Error in maxtree_BFilter: The aggregation "%s" does not exist.\n', aggregation);
    return;
end
if (~isfield(maxtree,infield_name))
    fprintf('Error in maxtree_BFilter: The field "%s" does not exist.\n', infield_name);
    return;
end

%% Loop on the tree nodes
if     (strcmp(filter,'Erosion') || ...
        strcmp(filter,'Dilation') || ...
        strcmp(filter,'Mean') || ...
        strcmp(filter,'Median') ...
        )
    lifo_queue = zeros(length(maxtree),2); % queue with nodes id and their distance for the filtered node
    for n=1:length(maxtree)
        %% Extract the values to be filtered
        n1_values = 1;
        values    = zeros(length(maxtree),1);
        values(n1_values) = maxtree(n).(infield_name);
        
        % Ancestors
        node = maxtree(n).Parent;
        for m=1:parameter
            if (isempty(node)) 
                break; 
            end
            n1_values = n1_values + 1;
            values(n1_values) = double(maxtree(node).(infield_name));
            node = maxtree(node).Parent;
        end
        
        % Descendants
        n2_values = 0;
        nb_values = 0;
        valuesB   = zeros(length(maxtree),1);
        lifo_stop = 0;
        if ((parameter ~= 0) && (length(maxtree(n).Children)~=0))
            % Introduce the children in the queue
            for i=1:length(maxtree(n).Children)
                lifo_stop = lifo_stop + 1;
                lifo_queue(lifo_stop,1) = maxtree(n).Children(i); 
                lifo_queue(lifo_stop,2) = 1; 
            end
  
            while(lifo_stop>0) 
                % Define the position of the node
                node = lifo_queue(lifo_stop,1); 
                dist = lifo_queue(lifo_stop,2); 
                lifo_stop = lifo_stop - 1;
                
                % Get the parent until node "n" if necessary
                if ((n2_values==0) && (dist~=1))
                    nodeTmp = maxtree(node).Parent;
                    while (nodeTmp ~= n)
                     	n2_values = n2_values + 1;
                        values(n1_values + n2_values) = double(maxtree(nodeTmp).(infield_name));
                        nodeTmp = maxtree(nodeTmp).Parent;
                    end
                end
                
                n2_values = n2_values + 1;
                values(n1_values + n2_values) = double(maxtree(node).(infield_name));
                             
                % Put the children in the queue 
                if (dist < parameter)
                    for i=length(maxtree(node).Children):-1:1
                        nodeC = maxtree(node).Children(i);
                        lifo_stop = lifo_stop + 1;
                        lifo_queue(lifo_stop,1) = nodeC;
                        lifo_queue(lifo_stop,2) = dist+1;
                    end
                end
                
                % Filter if we have enough samples or have reached a leaf 
                if ((isempty(maxtree(node).Children)) || (dist==parameter))
                    % Filter the values
                    valuesT   = values(1:n1_values+n2_values);
                    nb_values = nb_values + 1;
                    if (strcmp(filter,'Erosion'))
                        valuesB(nb_values) = min(valuesT);
                    elseif (strcmp(filter,'Dilation'))
                        valuesB(nb_values) = max(valuesT);
                    elseif (strcmp(filter,'Mean'))
                        valuesB(nb_values) = mean(valuesT);
                    elseif (strcmp(filter,'Median'))
                        valuesB(nb_values) = median(valuesT);
                    end
                    values(n1_values+1:n1_values+n2_values) = 0;
                    n2_values = 0;
                end
            end
            
            % Get the final value
            if (strcmp(aggregation,'Min'))
                maxtree_out(n).(outfield_name) = min(valuesB(1:nb_values));
            elseif (strcmp(aggregation,'Max'))
                maxtree_out(n).(outfield_name) = max(valuesB(1:nb_values));
            elseif (strcmp(aggregation,'Mean'))
                maxtree_out(n).(outfield_name) = mean(valuesB(1:nb_values));
            elseif (strcmp(aggregation,'Median'))
                maxtree_out(n).(outfield_name) = median(valuesB(1:nb_values));
            end

        end % loop on branches
        
        % In case the node to filter is a leaf
        if (nb_values == 0)
            % Get the final value
            if (strcmp(aggregation,'Min'))
                maxtree_out(n).(outfield_name) = min(values(1:n1_values));
            elseif (strcmp(aggregation,'Max'))
                maxtree_out(n).(outfield_name) = max(values(1:n1_values));
            elseif (strcmp(aggregation,'Mean'))
                maxtree_out(n).(outfield_name) = mean(values(1:n1_values));
            elseif (strcmp(aggregation,'Median'))
                maxtree_out(n).(outfield_name) = median(values(1:n1_values));
            end
        end
        
    end % loop on n

end