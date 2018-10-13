function  maxtree_out = maxtree_TFilter(maxtree, filter, parameter, infield_name, outfield_name)
%  MAXTREE_TFILTER computes a Tree Filter on attributes of Maxtrees. In
%  Tree filters, the neigborhood does not include the siblings of the node
%  to be filtered nor their descendants. 
%
%  maxtree_out = MAXTREE_TFILTER(maxtree, filter, parameter, infield_name, outfield_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     filter:               Type of filters, possible values are: 
%                           'Erosion','Dilation',
%                           'Opening','Closing',
%                           'Mean','Median'
%     parameter:            size of the neighborhood in the branch
%     infield_name:         Name of the field defining the input data
%     outfield_name:        Name of the field defining the output data
% 
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%     maxtree = MAXTREE_TFILTER(maxtree,'Opening',2,'GrayLevel','Out');
%
%  See also MAXTREE_GFILTER, MAXTREE_BFILTER
%
%  Author: Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
maxtree_out = maxtree;
if (~strcmp(filter,'Erosion')  && ...
    ~strcmp(filter,'Dilation') && ...  
    ~strcmp(filter,'Opening')  && ...
    ~strcmp(filter,'Closing')  && ...
    ~strcmp(filter,'Mean')  && ...
    ~strcmp(filter,'Median')...
    ) 
    fprintf('Error in maxtree_TFilter: The filter "%s" does not exist.\n', filter);
    return;
end
if (~isfield(maxtree,infield_name))
    fprintf('Error in maxtree_TFilter: The field "%s" does not exist.\n', infield_name);
    return;
end

%% Loop on the tree nodes
if     (strcmp(filter,'Erosion') || ...
        strcmp(filter,'Dilation') || ...
        strcmp(filter,'Mean') || ...
        strcmp(filter,'Median') ...
        )
    fifo_queue = zeros(length(maxtree),2); % queue with nodes id and their distance for the filtered node
    for n=1:length(maxtree)
        %% Extract the values to be filtered
        values = double(maxtree(n).(infield_name));
        % Ancestors
        node = maxtree(n).Parent;
        for m=1:parameter
            if (isempty(node)) break; end
            values = [values double(maxtree(node).(infield_name))];
            node = maxtree(node).Parent;
        end
        % Descendants
        fifo_start = 1; fifo_stop  = 1;
        if (parameter ~= 0)  % Introduce the children in the queue
            for i=1:length(maxtree(n).Children)
                fifo_queue(fifo_stop,1) = maxtree(n).Children(i); 
                fifo_queue(fifo_stop,2) = 1; 
                fifo_stop  = fifo_stop +1;
            end
            while(fifo_start < fifo_stop) 
                node = fifo_queue(fifo_start,1); 
                dist = fifo_queue(fifo_start,2);
                fifo_start = fifo_start +1;
                values = [values double(maxtree(node).(infield_name))];
                if (dist<parameter)
                    for i=1:length(maxtree(node).Children)
                        fifo_queue(fifo_stop,1) = maxtree(node).Children(i); 
                        fifo_queue(fifo_stop,2) = dist+1; 
                        fifo_stop  = fifo_stop +1;
                    end
                end
            end
        end

        %% Filter the values
        if (strcmp(filter,'Erosion'))
            maxtree_out(n).(outfield_name) = min(values);
        elseif (strcmp(filter,'Dilation'))
            maxtree_out(n).(outfield_name) = max(values);
        elseif (strcmp(filter,'Mean'))
            maxtree_out(n).(outfield_name) = mean(values);
        elseif (strcmp(filter,'Median'))
            maxtree_out(n).(outfield_name) = median(values);
        end

    end
elseif (strcmp(filter,'Opening'))
    maxtree     = maxtree_TFilter(maxtree,'Erosion' ,parameter,infield_name, 'tmp');
    maxtree_out = maxtree_TFilter(maxtree,'Dilation',parameter,'tmp', outfield_name);
elseif (strcmp(filter,'Closing'))
    maxtree     = maxtree_TFilter(maxtree,'Dilation' ,parameter,infield_name, 'tmp');
    maxtree_out = maxtree_TFilter(maxtree,'Erosion',parameter,'tmp', outfield_name);
end

end