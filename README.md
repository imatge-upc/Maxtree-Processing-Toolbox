
# Maxtree Processing Toolbox
The Maxtree Processing Toolbox is a Matlab toolbox allowing easy creation, processing and handling of Maxtree or Mintree representations. Maxtree and Mintree are image (or signal) representations created by structuring the connected components resulting from threshold decomposition [1], [2]. They are also known as Component Trees [3]. They provide a multiscale description of extremas of images and signals and are, among other things, one of the classical ways to build connected operators [4]. These trees can also be populated with attributes,  resulting in *Graph attribute signals* that can themselves be processed. The main goal of this toolbox is to provide functions for the creation of trees and graph attribute signals as well as their processing. Many tools included in this toolbox are further discussed in [5].

The main functions include: 

* **Create** Maxtree or Mintree from either images or Maxtrees & Mintrees,
* **Populate** the Max/Mintree with attributes to create graph attribute signals, 
* **Visualize** Max/Mintrees to analyze the tree structure or the attribute signal evolution,
* **Process** attribute signals defined on Max/Mintrees through:
  * **Pruning** (remove nodes of the tree)
  * **Filtering** (filter the attribute signal with mean, median and morphological operators)
  * **Morphological reconstruction**
  * **Regional extrema** identification and **labeling** of connected component on the tree structure. 
* **Restitution** of processed signals and attributes either on image or Max/Mintrees.

**A complete documentation can be found [here](Doc/Maxtree_Processing_Toolbox_Doc.pdf)**

This toolbox has been designed to explore processing ideas and perform experiments on small size images. It has not been particularly optimized to minimize the CPU load nor the memory consumption. Therefore, if large size images have to processed, a proper optimized C/C++ implementation of the functions presented here should be used.

----

**References**

 1. P. Salembier, A. Oliveras, and L. Garrido. Motion connected operators for image sequences, In VIII European Signal Processing Conference, EUSIPCO'96, pages 1083-1086, Trieste, Italy, September 1996.
 2. P. Salembier, A. Oliveras, and L. Garrido. Anti-extensive connected operators for image and sequence processing. IEEE Transactions on Image Processing, 7(4):555{570, April 1998.
 3. R. Jones. Component trees for image ltering and segmentation. In 1997 IEEE Workshop on Nonlinear Signal and Image Processing, Mackinac Island, USA, 1997.
 4. P. Salembier and M. H. F. Wilkinson, Connected operators: A review of region-based morphological image processing techniques, IEEE Signal Processing Magazine, vol. 6, pp. 136–157, 2009.
 5. P. Salembier, S. Liesegang and C. Lopez-Martinez, [Ship Detection in SAR Images based on Maxtree Representation and Graph Signal Processing](Doc/Ship_Detection_in_SAR_Images_based_on_Maxtree_Representation_and_Graph_Signal_Processing.pdf), IEEE Transactions on Geoscience and Remote Sensing, Accepted for publication, 2018.

**Contact:** 
Philippe Salembier (email: <philippe.salembier@upc.edu>), Sergi Liesegang.
