function  maxtree_out = maxtree_GFilter(maxtree, filter, parameter, infield_name, outfield_name)
%  MAXTREE_GFILTER computes a Graph Filter on attributes of Maxtrees. In
%  Graph filters, the neigborhood is defined on the basis of the graph
%  connectivity
%
%  maxtree_out = MAXTREE_GFILTER(maxtree, filter, parameter, infield_name, outfield_name);
% 
%  Input arguments:
%     maxtree:              Maxtree structure
%     filter:               Type of filters, possible values are: 
%                           'Erosion','Dilation',
%                           'Opening','Closing',
%                           'Mean','Median'
%     parameter:            Size of the neighborhood in the graph
%     infield_name:         Name of the field defining the input data
%     outfield_name:        Name of the field defining the output data
% 
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%     maxtree = MAXTREE_GFILTER(maxtree,'Opening',2,'GrayLevel','Out');
%
%  See also MAXTREE_TFILTER, MAXTREE_BFILTER
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
    fprintf('Error in maxtree_GFilter: The filter "%s" does not exist.\n', filter);
    return;
end
if (~isfield(maxtree,infield_name))
    fprintf('Error in maxtree_GFilter: The field "%s" does not exist.\n', infield_name);
    return;
end

%% Loop on the tree nodes
if     (strcmp(filter,'Erosion') || ...
        strcmp(filter,'Dilation') || ...
        strcmp(filter,'Mean') || ...
        strcmp(filter,'Median') ...
        )
    fifo_queue      = zeros(length(maxtree),2); % queue with nodes id and their distance for the filtered node
    for n=1:length(maxtree)

        % Define the neighborhood 
        fifo_start = 1; fifo_stop  = 1;
        fifo_queue(fifo_stop,1) = n; 
        fifo_queue(fifo_stop,2) = 0; 
        fifo_stop  = fifo_stop +1;
        maxtree(n).InQueue = 1;
        nvalues = 0;
        values = zeros(length(maxtree),1); 
        nodes  = zeros(length(maxtree),1);
        
        if (parameter ~= 0)              
            while(fifo_start < fifo_stop) 
                % Extract the node from the queue
                node = fifo_queue(fifo_start,1); 
                dist = fifo_queue(fifo_start,2);
                fifo_start = fifo_start +1;
                nvalues = nvalues + 1; 
                values(nvalues) = double(maxtree(node).(infield_name));
                nodes(nvalues)  = node;
                
                if (dist<parameter)
                    % Store the parent in the queue (if necessary)
                    nodeP = maxtree(node).Parent;
                    if (~isempty(nodeP) && isempty(maxtree(nodeP).InQueue))
                        fifo_queue(fifo_stop,1) = nodeP; 
                        fifo_queue(fifo_stop,2) = dist+1; 
                        fifo_stop  = fifo_stop +1;
                        maxtree(nodeP).InQueue = 1;
                    end
                    % Store the children in the queue (if necessary)
                    for i=1:length(maxtree(node).Children)
                        nodeC = maxtree(node).Children(i);
                        if (isempty(maxtree(nodeC).InQueue))
                            fifo_queue(fifo_stop,1) = nodeC; 
                            fifo_queue(fifo_stop,2) = dist+1; 
                            fifo_stop  = fifo_stop +1;
                            maxtree(nodeC).InQueue = 1;
                        end
                    end
                end
            end
        else
            nvalues = 1;
            values(nvalues) = double(maxtree(n).(infield_name));
        end

        %% Filter the values
        values_tmp = values(1:nvalues);
        if (strcmp(filter,'Erosion'))
            maxtree_out(n).(outfield_name) = min(values_tmp);
        elseif (strcmp(filter,'Dilation'))
            maxtree_out(n).(outfield_name) = max(values_tmp);
        elseif (strcmp(filter,'Mean'))
            maxtree_out(n).(outfield_name) = mean(values_tmp);
        elseif (strcmp(filter,'Median'))
            maxtree_out(n).(outfield_name) = median(values_tmp);
        end
        nodes_tmp = nodes(1:nvalues);
        for nn=1:length(nodes_tmp)
            maxtree(nodes(nn)).InQueue = [];
        end

    end
elseif (strcmp(filter,'Opening'))
    maxtree     = maxtree_GFilter(maxtree,'Erosion' ,parameter,infield_name, 'tmp');
    maxtree_out = maxtree_GFilter(maxtree,'Dilation',parameter,'tmp', outfield_name);
elseif (strcmp(filter,'Closing'))
    maxtree     = maxtree_GFilter(maxtree,'Dilation' ,parameter,infield_name, 'tmp');
    maxtree_out = maxtree_GFilter(maxtree,'Erosion',parameter,'tmp', outfield_name);
end

end