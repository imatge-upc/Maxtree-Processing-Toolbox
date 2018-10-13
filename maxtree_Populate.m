function  tree_out = maxtree_Populate(tree, attribute, varargin)
%  MAXTREE_POPULATE populates a mintree or a maxtree with an attribute.
%
%  tree_out = MAXTREE_POPULATE(tree, attribute);
%
%  Input arguments:
%     tree:                 Structure with the maxtree or mintree 
%     attribute:            Attribute to be computed. Available options:
%                           'Area', 'MeanGrayLevel', 'Ellipse',
%                           'AreaExtinction', 'ContrastExtinction'
%  Optional arguments:
%     image:                In the case of the 'Ellipse' attribute an 
%                           image has to be passed to the function.
%  Output argument:
%     tree_out:             Output maxtree with the resulting attribute
%                           field
%  EXAMPLE
%     tree_out = MAXTREE_POPULATE(tree, 'Area');
%     tree_out = MAXTREE_POPULATE(tree, 'Ellispe', image);
%
%  Author: Sergi Liesegang Maria, Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%% Check input parameters
tree_out = tree;
if (~strcmp(attribute,'Area')  && ...
    ~strcmp(attribute,'MeanGrayLevel') && ...
    ~strcmp(attribute,'Ellipse') && ...   
    ~strcmp(attribute,'AreaExtinction') && ...
    ~strcmp(attribute,'ContrastExtinction') ...
    ) 
    fprintf('Error in maxtree_Populate: The attribute "%s" is not implemented.\n', attribute);
    return;
end

%% Get additional parameter if necessary
if (strcmp(attribute,'Ellipse'))
    if (length(varargin)==0)
            fprintf('Error in maxtree_Populate: The image has to be passed to the function as third argument.\n');
    else
        image = varargin{1};
    end
end


%% Populate the tree_out
if strcmp(attribute,'Area')
    for node=length(tree):-1:1
        tree_out(node).(attribute) = int32(tree(node).NumberOfPixels);
        if (~isempty(tree(node).Children))
            tree_out(node).(attribute) = tree_out(node).(attribute) ...
                + sum([tree_out(tree(node).Children).(attribute)]);
        end
    end
elseif strcmp(attribute,'MeanGrayLevel')
    for node=length(tree):-1:1
        tree(node).Area     = tree(node).NumberOfPixels;
        tree(node).Integral = double(tree(node).GrayLevel) * double(tree(node).NumberOfPixels);
        if (~isempty(tree(node).Children))
            tree(node).Area     = tree(node).Area + sum([tree(tree(node).Children).Area]);
            tree(node).Integral = tree(node).Integral + sum([tree(tree(node).Children).Integral]);
        end
        tree_out(node).(attribute) = double(tree(node).Integral)/double(tree(node).Area);
    end
elseif strcmp(attribute,'Ellipse')
    % First compute the component pixels for each node
    for node=length(tree):-1:1
        tree(node).ComponentPixels = tree(node).Pixels;
        if (~isempty(tree(node).Children))
            tree(node).ComponentPixels = sort([tree(node).ComponentPixels tree(tree(node).Children).ComponentPixels]);
        end
    end
    % Define the binary max representing each component and define the
    % ellipse
    num_of_pixels   = size(image,1) * size(image,2);
    pixel_positions = int32(0:num_of_pixels-1);
    for node=1:length(tree)
        remove_pixels = pixel_positions;
        remove_pixels(tree(node).ComponentPixels+1) = [];
        vector_image  = pixel_positions;
        vector_image(remove_pixels+1) = -1;
        binary_mask = (vec2mat((pixel_positions == vector_image),size(image,1)))';
        ellipse = regionprops(binary_mask,'MajorAxisLength',...
            'MinorAxisLength','Orientation','Eccentricity','BoundingBox');
        tree_out(node).Ellipse = ellipse;
    end
elseif strcmp(attribute,'AreaExtinction')
    fifo_queue = zeros(2*length(tree),1);
    fifo_start = 1; fifo_stop  = 1;
    tmp =  num2cell(zeros(1,length(tree)));
    [tree(:).Processed]          = tmp{:};
    [tree_out(:).AreaExtinction] = tmp{:};

    % Put all leaves in the queue
    for node=1:length(tree)
        if (length(tree(node).Children) == 0)
            fifo_queue(fifo_stop) = node; 
            fifo_stop  = fifo_stop +1;
        end
    end
    % Propagation from leaves to root
    while(fifo_start < fifo_stop) 
        node = fifo_queue(fifo_start); 
        fifo_start = fifo_start +1;

        % If there is just one branch going through the node
        if (length(tree(node).Children) <= 1)
            tree(node).Area = int32(tree(node).NumberOfPixels);
            if (~isempty(tree(node).Children)) % Then there is just one child
                tree(node).Area = tree(node).Area + tree(tree(node).Children).Area;
            end
            nodeP = tree(node).Parent;
            if (~isempty(nodeP))
                fifo_queue(fifo_stop) = nodeP; 
                fifo_stop  = fifo_stop +1;
                tree(node).Processed = 1;
            end
        else % If there are various branches starting at the node
            % Have all children been processed? 
            AllProcessedChildren = prod([tree(tree(node).Children).Processed]);
            
            if ((tree(node).Processed==0) && (AllProcessedChildren==1))
                tree(node).Area = int32(tree(node).NumberOfPixels);

                % Define the extinction values of children (minus one)
                maxArea = max([tree(tree(node).Children).Area]);
                for i=1:length(tree(node).Children)
                    nodeC = tree(node).Children(i);
                    if (tree(nodeC).Area == maxArea)
                        tree(node).Area = tree(node).Area + tree(nodeC).Area;
                    else
                        tree_out(nodeC).AreaExtinction = tree(nodeC).Area;
                    end
                end
                
                % Continue the propagation
                tree(node).Processed = 1;
                nodeP = tree(node).Parent;
                if (~isempty(nodeP))
                    fifo_queue(fifo_stop) = nodeP; 
                    fifo_stop  = fifo_stop +1;

                end
            end
        end
    end
    % Define the extinction value for the root node
    tree_out(1).AreaExtinction = tree(1).Area;
    
    % Propagate the extinction values to all nodes (root to leaves)
    for node=2:length(tree_out)
        if (tree_out(node).AreaExtinction == 0)
            tree_out(node).AreaExtinction = tree_out(tree_out(node).Parent).AreaExtinction;
        end
    end
elseif strcmp(attribute,'ContrastExtinction')
    fifo_queue = zeros(2*length(tree),1);
    fifo_start = 1; fifo_stop  = 1;
    tmp =  num2cell(zeros(1,length(tree)));
    [tree(:).Processed]              = tmp{:};
    [tree_out(:).ContrastExtinction] = tmp{:};

    % Put all leaves in the queue
    for node=1:length(tree)
        if (length(tree(node).Children) == 0)
            fifo_queue(fifo_stop) = node; 
            fifo_stop  = fifo_stop +1;
        end
    end
    % Propagation from leaves to root
    while(fifo_start < fifo_stop) 
        node = fifo_queue(fifo_start); 
        fifo_start = fifo_start +1;

        % If there is just one branch going through the node
        if (length(tree(node).Children) <= 1)
            tree(node).Contrast = 0;
            if (~isempty(tree(node).Children)) % Then there is just one child
                tree(node).Contrast = tree(tree(node).Children).GrayLevel ...
                    - tree(node).GrayLevel + tree(tree(node).Children).Contrast;
            end
            nodeP = tree(node).Parent;
            if (~isempty(nodeP))
                fifo_queue(fifo_stop) = nodeP; 
                fifo_stop  = fifo_stop +1;
                tree(node).Processed = 1;
            end
        else % If there are various branches starting at the node
            % Have all children been processed? 
            AllProcessedChildren = prod([tree(tree(node).Children).Processed]);
            
            if ((tree(node).Processed==0) && (AllProcessedChildren==1))
                tree(node).Contrast = 0;

                % Define the extinction values of children (minus one)
                maxContrast = 0;
                for i=1:length(tree(node).Children)
                    nodeC = tree(node).Children(i);
                    CurrentContrast = tree(nodeC).GrayLevel ...
                        - tree(node).GrayLevel + tree(nodeC).Contrast;
                    maxContrast = max(maxContrast, CurrentContrast);
                end
                for i=1:length(tree(node).Children)
                    nodeC = tree(node).Children(i);
                    CurrentContrast = tree(nodeC).GrayLevel ...
                        - tree(node).GrayLevel + tree(nodeC).Contrast;
                    if (CurrentContrast == maxContrast)
                        tree(node).Contrast = CurrentContrast;
                    else
                        tree_out(nodeC).ContrastExtinction = tree(nodeC).GrayLevel ...
                            - tree(node).GrayLevel + tree(nodeC).Contrast;
                    end
                end
                
                % Continue the propagation
                tree(node).Processed = 1;
                nodeP = tree(node).Parent;
                if (~isempty(nodeP))
                    fifo_queue(fifo_stop) = nodeP; 
                    fifo_stop  = fifo_stop +1;

                end
            end
        end
    end
    % Define the extinction value for the root node
    tree_out(1).ContrastExtinction = tree(1).Contrast;
    
    % Propagate the extinction values to all nodes (root to leaves)
    for node=2:length(tree_out)
        if (tree_out(node).ContrastExtinction == 0)
            tree_out(node).ContrastExtinction = tree_out(tree_out(node).Parent).ContrastExtinction;
        end
    end
end

end