% /* ==========================================================
%  * UNIVERSITAT POLITECNICA DE CATALUNYA, Spain
%  * Image processing group: https://imatge.upc.edu/web/
%  * ==========================================================
%  *
%  * mintree_of_image.m 
%  *
%  * Compute the mintree representation of a matlab int16 image 
%  * and outputs matlab mintree as a maxtree structure, which is a 
%  * struct array with the following fields: 
%  * Node:            ID on the node in the tree
%  * GrayLevel:       Gray level of the pixels stored in the node
%  * Parent:          ID of the Parent nodes
%  * Children:        ID(s) of the Child node
%  * NumberOfPixels   Number of pixels stored inthe node
%  * Pixels           Offest from the base defining the pixel position
%  *
%  * The calling syntax is:
%  *      mintree = mintree_of_image(image, connectivity)
%  *
%  *      image should be an int16 array whose values are in [0,32000]
%  *      connectivity should be either 4 or 8 
%  *
%  * Authors of the maxtree related functions : 
%  *    Luis Garrido Ostermann, Philippe Salembier
%  *
%  * References on Maxtree representation: 
%  * [1]  P. Salembier, Oliveras, A., and Garrido, L., Motion connected 
%  *      operators for image sequences, in VIII European Signal Processing 
%  *      Conference, EUSIPCO'96, Trieste, Italy, 1996, pp.1083-1086.
%  * [2]  P. Salembier, Oliveras, A., and Garrido, L., Antiextensive 
%  *      connected operators for image and sequence processing, 
%  *      IEEE Transactions on Image Processing, vol.7, pp.555-570,1998
%  *
%  *========================================================*/

function  mintree = mintree_of_image(ima, connectivity)

% Compute the mintree by duality
maxval = 32000;
mintree = maxtree_of_image(maxval-ima,connectivity);

for n=1:length(mintree)
    mintree(n).GrayLevel = maxval - mintree(n).GrayLevel;
end

end
