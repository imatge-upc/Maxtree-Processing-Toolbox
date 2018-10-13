%  Examples of function use of the Maxtree Processing Toolbox
%
%  Example 1: Simple maxtree construction, area filtering and image restitution 
%  Example 2: Maxtree construction, population with a random feature,
%             graph, tree and branch filtering of the feature, trees 
%             visualization with graphviz
%  Example 3: Extract the regional maxima of a feature on a tree and label them, 
%             tree visualization with graphviz
%  Example 4: Filter by reconstruction (Branch and Tree) applied on maxtree 
%             (markers are the absolute maxima of the feature), trees 
%             visualization with graphviz
%  Example 5: Compute a maxtree of a maxtree
%  Example 6: Compute a maxtree of a maxtree, AreaExtinction Filter the second 
%             maxtree to remove small maxima, Restore the filtered features in 
%             the first maxtree
%  Example 7: Visualization example: Construct a maxtree, visualization 
%             with maxtree_to_gv (Graphviz), maxtree_branches_display and maxtree_plot 
%
%  Authors:   Philippe Salembier 
%  Copyright 2016, UPC, Image processing group, https://imatge.upc.edu

%%
clear;
close all;
Example = 1;

%%
if (Example==1)
    % Creation of the image and of the maxtree
    ima = int16(   [1 1 1 1 1 1;
                    1 1 3 1 1 2;
                    4 1 2 3 1 3;
                    2 6 1 1 1 4;
                    6 5 3 1 1 5])
    tree = maxtree_of_image(ima,4);
    
    % Tree population with area feature and pruning
    tree = maxtree_Populate(tree,'Area');
    tree = maxtree_Prune(tree,'Area',3,'Direct');
    
    % Restitution of the image after area opening
    ima_out = image_Restitution(tree,'GrayLevel',size(ima,1),size(ima,2));
    ima_out
end

%%
if (Example==2)
    % Filter definition
    Filter = 'Median';
    Parameter = 2;
    
    % Creation of the image and maxtree
    ima = int16(   [1 1 1 1 1;
                    1 1 3 1 1;
                    4 1 2 3 1;
                    2 6 1 1 1;
                    6 5 3 1 1]);

    tree = maxtree_of_image(ima,4);
    [branches, tree] = maxtree_to_branches(tree);
    
    % Tree population with a random feature
    rng(4);
    tmp =  num2cell(int32(8*rand(1,length(tree))));
    [tree(:).Feature] = tmp{:};
    maxtree_to_gv('maxtree.gv','dot', tree,'Attribute','Feature','Bounds',[0 8])    

    % Graph filtering 
    tree = maxtree_GFilter(tree, Filter, Parameter, 'Feature', 'GFFeature');
    maxtree_to_gv('maxtree_GFilter.gv','dot', tree,'Attribute','GFFeature','Bounds',[0 8]);
    % Branch filtering 
    tree = maxtree_TFilter(tree, Filter, Parameter, 'Feature', 'TFFeature');
    maxtree_to_gv('maxtree_TFilter.gv','dot', tree,'Attribute','TFFeature','Bounds',[0 8]);
    % Tree filtering 
    tree = maxtree_BFilter(tree, Filter, Parameter, Filter, 'Feature', 'BFFeature');
    maxtree_to_gv('maxtree_BFilter.gv','dot', tree,'Attribute','BFFeature','Bounds',[0 8]);
    % Tree filtering (Alternative implementation)
    % tree = maxtree_BFilter2(tree, branches, Filter, Parameter, Filter, 'Feature', 'BFFeature2');
    % maxtree_to_gv('maxtree_BFilter2.gv','dot', tree,'Attribute','BFFeature2','Bounds',[0 8]);
end

%%
if (Example==3)
    % Creation of the image and maxtree
    ima = int16(   [1 1 1 1 1;
                    1 1 3 1 1;
                    4 1 2 3 1;
                    2 6 1 1 1;
                    6 5 3 1 1]);
    tree = maxtree_of_image(ima,4);
    
    % Tree population with a random feature
    rng(1);
    tmp =  num2cell(int32(8*rand(1,length(tree))));
    [tree(:).Feature] = tmp{:};
    maxtree_to_gv('maxtree.gv','dot', tree,'Attribute','Feature','Bounds',[0 8])    

    % Extract the Feature maxima and visualize
    tree = maxtree_RegExtrema(tree, 'Max', 'Feature', 'FMax');
    % Label the individual maximas
    tree = maxtree_Label(tree, 'FMax', 'FMaxLabel');

    maxtree_to_gv('maxtree_RMax.gv','dot',tree,'Attribute','FMaxLabel','DisplayMode','Color');
  
end

%%
if (Example==4)
    % Creation of the image and maxtree
    ima = int16(   [1 1 1 1 1;
                    1 1 3 1 1;
                    4 1 2 3 1;
                    2 6 1 1 1;
                    6 5 3 1 1]);
    tree = maxtree_of_image(ima,4);
    
    % Tree population with a random feature
    rng(14);
    tmp =  num2cell(int32(8*rand(1,length(tree))));
    [tree(:).Feature] = tmp{:};
    maxtree_to_gv('maxtree.gv','dot', tree,'Attribute','Feature','Bounds',[0 8])    

    % Create the marker (absolute max of the image) 
    marker =[tree(:).Feature];
    maxval = max(marker);
    marker(marker~=maxval)=0;
    marker = num2cell(marker);
    [tree(:).Marker] = marker{:};
   
    % Reconstruct (both Graph and Branch reconstructions) 
    tree = maxtree_TReconstruct(tree, 'Leaves_to_Root', 'Feature', 'Marker', 'Feature_TRec');
    tree = maxtree_GReconstruct(tree, 'Feature', 'Marker', 'Feature_GRec');

    maxtree_to_gv('maxtree_TRec.gv','dot', tree,'Attribute','Feature_TRec','Bounds',[0 8]);
    maxtree_to_gv('maxtree_GRec.gv','dot', tree,'Attribute','Feature_GRec','Bounds',[0 8]);

end

%%
if (Example==5)
    % Creation of the image and maxtree
    ima = int16(   [1 1 1 1 1;
                    1 1 3 1 1;
                    4 1 2 3 1;
                    2 6 1 1 1;
                    6 5 3 1 1]);
    tree1 = maxtree_of_image(ima,4);
    
    % Tree population with a random feature
    rng(5);
    tmp =  num2cell(int32(8*rand(1,length(tree1))));
    [tree1(:).Feature] = tmp{:};
    maxtree_to_gv('maxtree1.gv','dot', tree1,'Attribute','Feature','Bounds',[0 8]);    

    % Compute the maxtree of the maxtree
    tree2 = maxtree_of_maxmintree(tree1, 'Feature');
    maxtree_to_gv('maxtree2.gv','dot', tree2,'Attribute','GrayLevel','Bounds',[0 8]);

end

%%
if (Example==6)
    rng(1);
    ima = int16(8*rand(10,10));
    tree1 = maxtree_of_image(ima,4);
    
    % Tree population with a random feature
    rng(1);
    tmp =  num2cell(int32(7*rand(1,length(tree1))));
    [tree1(:).Feature] = tmp{:};
    
    maxtree_to_gv('maxtree1.gv','dot', tree1,'Attribute','Feature','Bounds',[0 4]);

    % Compute the maxtree of the maxtree
    tree2 = maxtree_of_maxmintree(tree1, 'Feature');
    maxtree_to_gv('maxtree2.gv','dot', tree2,'Attribute','GrayLevel','Bounds',[0 4]);

    % Tree population with area extinction feature (note that "area
    % extinction" and "area" are two different features) and pruning
    tree2 = maxtree_Populate(tree2,'AreaExtinction');
    tree2 = maxtree_Prune(tree2,'AreaExtinction',3,'Direct');
    maxtree_to_gv('maxtree2_filt.gv','dot', tree2,'Attribute','GrayLevel','Bounds',[0 4]);

    % Restitution of the inital maxtree after area opening
    tree3 = tree_Restitution(tree2,'GrayLevel',tree1,'GrayLevel_Processed');
    maxtree_to_gv('maxtree1_filt.gv','dot', tree3,'Attribute','GrayLevel_Processed','Bounds',[0 4]);

end

%%
if (Example==7)
    % Creation of the image and maxtree
    rng(1);
    ima = int16(8*rand(20,20));
    imshow(ima,[]);
    tree = maxtree_of_image(ima,4);
    [branches, tree] = maxtree_to_branches(tree);

    maxtree_to_gv('maxtree.gv','dot', tree);
    BranchGrayLevel = maxtree_branches_display(tree,branches,'GrayLevel');
    figure; imshow(BranchGrayLevel,[]);
    
    figure;
    maxtree_plot(tree, 100, 'GrayLevel', '-+b');

end