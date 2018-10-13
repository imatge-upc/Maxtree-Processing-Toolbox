function  maxtree_out = maxtree_BFilter2(maxtree, branches, filter, parameter, aggregation, infield_name, outfield_name)
%  MAXTREE_BFILTER2 computes a Branch Filter on attributes of Maxtrees. In
%  Branch filters, each branch is individually filtered (estimation step)
%  and then the resulting estimated values are aggregated. This function is
%  an alternative implementation of maxtree_BFilter. It is faster in
%  particular for large values of 'parameter' (Note that the padding
%  strategy is slightly different in both functions. So the results are
%  slightly different). 
% 
%  maxtree_out = MAXTREE_BFILTER2(maxtree, branches, filter, parameter,... 
%                aggregation, infield_name, outfield_name);
%
%  Input arguments:
%     maxtree:              Maxtree structure
%     branches:             Branch representation of the maxtree computed 
%                           with maxtree_to_branches
%     filter:               Type of (estimation) filter, possible values are: 
%                           'Erosion','Dilation', 'Mean','Median'
%     parameter:            size of the neighborhood in the branch
%     aggregation:          Type of aggregation, possible values are:
%                           'Mean','Median', 'Max', 'Min''
%     infield_name:         Name of the field defining the input data
%     outfield_name:        Name of the field defining the output data
%
%  Output argument:
%     maxtree_out:          Output maxtree with the resulting outfield_name
%                           field
%  EXAMPLE 
%     maxtree = MAXTREE_BFILTER2(maxtree,branches,'Opening',2,'Median','GrayLevel','Out');
%
%  See also MAXTREE_TO_BRANCHES, MAXTREE_GFILTER, MAXTREE_TFILTER, MAXTREE_BFILTER
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
    fprintf('Error in maxtree_BFilter: The (estimation) filter "%s" does not exist.\n', filter);
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

%% Estimation step
% Create the image to filter 
h = size(branches,1);
w = size(branches,2);
image1 = zeros(h+2*parameter,w);
image1(:,:) = inf;
for i=1:h
    for j=1:w
        if(branches(i,j)>0)
            image1(i+parameter,j) = maxtree(branches(i,j)).(infield_name);
        end
    end
end

% Mirroring below the root
for i=h+parameter+1:1:h+2*parameter
    for j=1:w
        if (image1(2*h+2*parameter+1-i,j)~=Inf)
            image1(i,j) = image1(2*h+2*parameter+1-i,j);
        else
            image1(i,j) = image1(i-1,j);
        end
    end
end
% Mirroring above the leaves
for j=1:w
    i=h+parameter;
    while ((i>=1)&&(image1(i,j)~=Inf))
        i=i-1;
    end
    i0=i;
    while (i>=1)
        if ((2*i0-i+1<=size(image1,1)) && (image1(2*i0-i+1,j)~=Inf))
            image1(i,j) = image1(2*i0-i+1,j);
        else
            image1(i,j) = image1(i+1,j);
        end
        i=i-1;
    end
end

% Estimation filtering
h_filter = ones(2*parameter+1,1);
if (strcmp(filter,'Erosion'))
    image1 = ordfilt2(image1,1,h_filter);
elseif (strcmp(filter,'Dilation'))
    image1 = ordfilt2(image1,2*parameter+1,h_filter);
elseif (strcmp(filter,'Median'))
    image1 = ordfilt2(image1,parameter+1,h_filter);
elseif (strcmp(filter,'Mean'))
    image1 = conv2(image1,h_filter/(2*parameter+1),same,'same');
end
image2 = image1(parameter+1:parameter+h,1:w);

%% Aggregation step
for n=1:length(maxtree)
    [row,col] = find(branches==n);
    values    = zeros(length(row),1);
    ivalue    =  1;
    node_end0 = -1;
    for i=1:length(row)
        node_end1 = branches(max(1,row(i)-parameter),col(i));
        if ((node_end1~=node_end0)||(node_end1==0))
            values(ivalue) = image2(row(i),col(i));
            ivalue = ivalue+1;
        end
        node_end0 = node_end1;
    end
    values = values(1:ivalue-1);
    % Get the final value
    if (strcmp(aggregation,'Min'))
        maxtree_out(n).(outfield_name) = min(values);
    elseif (strcmp(aggregation,'Max'))
        maxtree_out(n).(outfield_name) = max(values);
    elseif (strcmp(aggregation,'Mean'))
        maxtree_out(n).(outfield_name) = mean(values);
    elseif (strcmp(aggregation,'Median'))
        maxtree_out(n).(outfield_name) = median(values);
    end
end

end