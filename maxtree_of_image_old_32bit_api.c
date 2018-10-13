/* ==========================================================
 * UNIVERSITAT POLITECNICA DE CATALUNYA, Spain
 * Image processing group: https://imatge.upc.edu/web/
 * ==========================================================
 *
 * maxtree_of_image.c - MATLAB mex Interface
 *
 * Compute the maxtree representation of a matlab int16 image 
 * and outputs matlab maxtree structure, which is a struct array 
 * with the following fields: 
 * Node:            ID on the node in the tree
 * GrayLevel:       Gray level of the pixels stored in the node
 * Parent:          ID of the Parent nodes
 * Children:        ID(s) of the Child node
 * NumberOfPixels   Number of pixels stored inthe node
 * Pixels           Offest from the base defining the pixel position
 *
 * The calling syntax is:
 *		maxtree = maxtree_of_image(image, connectivity)
 *
 *      image should be an int16 array whose values are in [0,32000]
 *      connectivity should be either 4 or 8 
 *
 * Authors of the maxtree related functions : 
 *    Luis Garrido Ostermann, Philippe Salembier
 *
 * Part of this code was developped in the context of the 
 * European RACE Morpheco project:
 *    UNIVERSITAT POLITECNICA DE CATALUNYA, Spain
 *    PARIS SCHOOL OF MINES, ARMINES, France
 *    ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland
 *    IBERMATICA, Spain
 *    LABORATOIRE D'ELECTRONIQUE PHILIPS, France
 *    ZENON S.A., Greece
 *
 * References on Maxtree representation: 
 * [1]  P. Salembier, Oliveras, A., and Garrido, L., Motion connected 
 *      operators for image sequences, in VIII European Signal Processing 
 *      Conference, EUSIPCO'96, Trieste, Italy, 1996, pp.1083-1086.
 * [2]  P. Salembier, Oliveras, A., and Garrido, L., Antiextensive 
 *      connected operators for image and sequence processing, 
 *      IEEE Transactions on Image Processing, vol.7, pp.555-570,1998
 *
 *========================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "mex.h"
#include "maxtree_of_image.h"
#include "matrix.h"


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize          connectivity, num_nodes;/* Connectivity and number of nodes in the tree */
    short int       *pr, *pi;               /* Pointers to image data */
    mwSize          m, n;                   /* Image size */
    struct image    *ima;                   /* Image structure */
    short int       max_ima, min_ima;       /* bounds of the image gray levels */
    struct max_tree *tree;                  /* Maxtree C structure */
    char            **fnames;               /* field names of the matlab structure */ 
    mxArray         *fout;
    mwSize          i, j ,k, ID, tmp1, tmp2, tmp3, *tmp0;
    
    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:nlhs","One output required.");
    }
    
    /* make sure the first input argument is type Int16 */
    if( !mxIsInt16(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:notInt16","Input image must be type int16.");
    }
    
    /* make sure the second input argument is scalar */
    if( mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:notScalar","Connectivity should be a scalar");
    }
    /* get the value of the connectivity */
    connectivity = mxGetScalar(prhs[1]);
    if( (connectivity==8 || connectivity==4) !=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:Connectivity","Connectivity should be either 4 or 8");
    }

    /* Get the size of the image and transfer in a struct image */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    ima = (struct image *) malloc(sizeof(struct image));
    ima->version=2; ima->x=m; ima->y=n; ima->grid='s';
    ima->f = (short int *) calloc(ima->x*ima->y,(long) sizeof(short int));
    pi = ima->f;
    pr = (short int *)mxGetData(prhs[0]);
    max_ima = SMLVAL;
    min_ima = BIGVAL;
    for (k=0; k<m*n; k++){
        *pi = *pr; 
        if (*pi>max_ima) {max_ima = *pi;}
        if (*pi<min_ima) {min_ima = *pi;}
        pi++, pr++;}
    /* Check if the image bounds are correct [0,32000] */
    if( (max_ima>32000 || min_ima<0)==1 ) {
        mexErrMsgIdAndTxt("MyToolbox:maxtree_creation:ImageBounds","The image values should be in [0,32000]");
    }

    /* Call the maxtree_creation */
    tree = (struct max_tree *) malloc(sizeof(struct max_tree));
    create_max_tree(tree, ima, connectivity);
    
    /* Set the nodeID. This ID will determine the node position in the matlab structure */
    k = 1; /* nodeID initialization */
    for (i=0; i<GREY_LEVELS; i++){
        if (tree->level[i].nnodes!=0) {
            for (n=0; n<tree->level[i].nnodes; n++) {
                tree->level[i].node[n].nodeID = k;
                k++;
            }
        }
    }

    /* Compute the number of nodes */
    num_nodes = 0;
    for (i=0; i<GREY_LEVELS; i++){
        num_nodes = num_nodes + tree->level[i].nnodes;
        }
    // mexPrintf("Number on nodes: %d  \n",num_nodes);
    
    /* Create the output Matlab max_tree structure */
    fnames = (char**)malloc(6*sizeof(char*));
    fnames[0] = "Node";
    fnames[1] = "GrayLevel";
    fnames[2] = "Parent";
    fnames[3] = "Children";
    fnames[4] = "NumberOfPixels";
    fnames[5] = "Pixels";
    plhs[0] = mxCreateStructMatrix((mwSize)num_nodes,1,6,(const char **)fnames);
    
    /* Transfer the C structure to the Matlab structure */
    for (i=0; i<GREY_LEVELS; i++){
        if (tree->level[i].nnodes!=0) {
            for (n=0; n<tree->level[i].nnodes; n++) {
                /* Write the nodeID */ 
                ID   = (mwSize)tree->level[i].node[n].nodeID;
                fout = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
                tmp0 = (mwSize *) mxGetPr(fout); tmp0[0] = ID;
                mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 0, fout);

                /* Write the Gray level */ 
                // fout = mxCreateDoubleScalar((double)i);
                fout = mxCreateNumericMatrix(1, 1, mxINT16_CLASS, mxREAL);
                tmp0 = (mwSize *) mxGetPr(fout); tmp0[0] = i;
                mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 1, fout);

                /* Write the Parent ID */
                // mexPrintf("             Parent (level,node): %d, %d  \n",tree->level[i].node[n].father.level,tree->level[i].node[n].father.node);
                tmp1 = tree->level[i].node[n].father.level;
                tmp2 = tree->level[i].node[n].father.node;
                fout = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
                if (tmp1!=-1){
                    tmp3 = tree->level[tmp1].node[tmp2].nodeID;
                    tmp0 = (mwSize *) mxGetPr(fout); tmp0[0] = tmp3;
                    mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 2, fout);
                }

                /* Write the Children IDs */
                // mexPrintf("             Number of children: %d \n",tree->level[i].node[n].NChildren);
                tmp1    = tree->level[i].node[n].NChildren;   if (tmp1<1) {tmp1=1;}
                fout    = mxCreateNumericMatrix(1, tmp1, mxINT32_CLASS, mxREAL);
                tmp0    = (mwSize *) mxGetPr(fout);                  
                if (tree->level[i].node[n].NChildren>0){
                    for (m=0; m<tree->level[i].node[n].NChildren; m++){
                        // mexPrintf("             Children: %d, (level,node): %d, %d \n",m,tree->level[i].node[n].Child[m].level, tree->level[i].node[n].Child[m].node);
                        tmp1 = tree->level[i].node[n].Child[m].level;
                        tmp2 = tree->level[i].node[n].Child[m].node;
                        tmp3 = tree->level[tmp1].node[tmp2].nodeID;
                        tmp0[m] = tmp3;
                    }
                mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 3, fout);
                }
       
                /* Write the Number of Pixels */ 
                // mexPrintf("    node: %d, Number of pixels: %d  \n",n,(int)tree->level[i].node[n].NPixels);
                tmp1 = (mwSize)tree->level[i].node[n].NPixels;
                // fout = mxCreateDoubleScalar((double)tmp1);
                fout = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
                tmp0 = (mwSize *) mxGetPr(fout); tmp0[0] = tmp1;
                mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 4, fout);

                /* Write the Pixel offsets */ 
                // for (m=0; m<tree->level[i].node[n].NPixels; m++){
                //    mexPrintf("             Pixel: %d, %d coord:(%d,%d)\n",m,(int)tree->level[i].node[n].Pixel[m],
                //           (int)PIX_COOR_X(ima,tree->level[i].node[n].Pixel[m]),
                //           (int)PIX_COOR_Y(ima,tree->level[i].node[n].Pixel[m]));}
                tmp1    = tree->level[i].node[n].NPixels;   
                fout    = mxCreateNumericMatrix(1, tmp1, mxINT32_CLASS, mxREAL);
                tmp0    = (mwSize *) mxGetPr(fout);                  
                for (m=0; m<tree->level[i].node[n].NPixels; m++){
                    tmp3 = tree->level[i].node[n].Pixel[m];
                    tmp0[m] = tmp3;
                }
                mxSetFieldByNumber(plhs[0], (mwIndex)(ID-1), 5, fout);
                
            }
        }
      }
}


/*****************************************************************************/
/*
 *
 *  This section includes functions to create and destroy Max Trees.
 *  The Max Tree structures are defined in *max_tree.h*. 
 *  The nodes are allocated in a matrix-like structure.
 *  struct max_tree tree;           // Suppose that the tree is now created
 *  tree->level[100].nnodes          --> number of nodes associated to the 
 *                                       gray-level 100.
 *  tree->level[100].node[2]         --> access 3rd node of gray-level 100.
 *  tree->level[100].node[2].father  --> coordinates (gray-level, nuber-node)
 *                                       of its father
 *  tree->level[100].node[2].NPixels --> number of pixels associated to the
 *                                       node (100,2).
 ******************************************************************************/

/******************************************************************************/
/*
 *  create_max_tree
 */
/*    Create the Max Tree representation of image "img" with
 *    connectivity "connectivity" (4 or 8). Image gray levels should
 *    be between 0 and 32000.
 *
 */

#define NOT_ANALYSED         -1
#define IN_QUEUE             -2
#define BORDER               BIGVAL

void create_max_tree (struct max_tree *tree, struct image *img, int connectivity)
{
    struct image *map;
    struct hqueue *list_to_propagate;
    struct max_node *p_node;
    long j, i, end;
    datim *min_value, *node_at_level;
    int direction[8];
    
    
    /* Here begins code of create_max_tree */
    
    if ((connectivity!=8)&&(connectivity!=4))
        runerr("ERROR(Create_tree) : Connectivity must be 4 or 8",1);
    
    /* Init the tree */
    
    tree->image_sizx = img->x;
    tree->image_sizy = img->y;
    
    tree->level = (struct array_max_nodes *)
    emalloc((sizs) (sizeof(struct array_max_nodes) * GREY_LEVELS));
    
    for(j = 0; j<GREY_LEVELS; j++)
    {
        tree->level[j].nnodes = -1;
        tree->level[j].node = 0x0;
    }
    
    /* Init the Hqueue */
    
    list_to_propagate = new_fifo_hqueue(GREY_LEVELS - 1, 0);
    
    /* Init map (aux image) */
    
    map = (struct image *) emalloc((sizs) sizeof(struct image));
    timmalloc(img, map);
    imera(map, NOT_ANALYSED);
    
    putbrd(img, BORDER);
    putbrd(map, BORDER);
    
    /* pointers have been reallocated */
    
    direction[0] = 1; direction[1] = (int) img->x;
    direction[2] = -1; direction[3] = (int) -img->x;
    direction[4] = 1 + (int) img->x; direction[5] = (int) img->x - 1;
    direction[6] = -1 - (int) img->x; direction[7] = 1 - (int) img->x;
    
    /* we calculate *minvalue */
    
    min_value = img->f;
    end = (img->x)*(img->y);
    for(j = 0; j < end; j++) min_value=(*min_value<*(img->f+j))?min_value:(img->f+j);
    
    /* Init node_at_level */
    
    node_at_level = (datim *) ecalloc((sizs) GREY_LEVELS, (sizs) sizeof(datim));
    node_at_level[*min_value] = TRUE;
    
    /* Init for propagation */
    
    putvqueue(min_value, list_to_propagate);
    *(map->f + (min_value - img->f)) = IN_QUEUE;
    
    /* Let's propagate */
    
    flood(*min_value, list_to_propagate, tree, map, img, connectivity, direction,
          node_at_level, p_node);
    /* Ok. Let's reallocate ... */
    
    for(j = 0; j<GREY_LEVELS; j++)
    {
        if (++tree->level[j].nnodes != 0)
            tree->level[j].node = (struct max_node *)
            realloc(tree->level[j].node,
                    sizeof(struct max_node) * (tree->level[j].nnodes));
    }
    
    /* Compute children */
    
    compute_children(tree);
    
    /* ... and deallocate */
    
    free_hqueue(list_to_propagate);
    free((char *) list_to_propagate);
    
    rmbrd(img);
    rmbrd(map);
    
    /* Scan image in order to count NPixels of each node */
    
    end = (img->x) * (img->y);
    for(j = 0; j<end; j++)
    {
        p_node = tree->level[*(img->f+j)].node + (*(map->f + j));
        p_node->NPixels++;
    }
    
    /* Allocate memory for pixels */
    
    for(j = 0; j<GREY_LEVELS; j++)
    {
        end = tree->level[j].nnodes;
        for(i = 0; i<end; i++)
        {
            p_node = tree->level[j].node + i; /* p_node = &(tree->level[j].node[i]) */
            p_node->Pixel = (long *) emalloc((sizs) (sizeof(long) * p_node->NPixels));
            p_node->NPixels = 0; /* We'll use it as a counter */
        }
    }
    
    /* Allocate pixels (offset from base address) */
    
    end = (img->x) * (img->y);
    for(j = 0; j<end; j++) 
    {
        p_node = tree->level[*(img->f+j)].node + (*(map->f+j));
        p_node->Pixel[p_node->NPixels++] = j;
    }
    
    /* Deallocate map */    
    
    imfree(map);
    free((char *) map);
    
    free((char *) node_at_level);
}

/* Floding proces */

void flood(datim lev, struct hqueue *list_to_propagate, struct max_tree *tree,
     struct image *map, struct image *img, int connectivity, int *direction,
     datim *node_at_level, struct max_node *p_node)
 {
   datim *neighbour_pixel, *pixel, l, index, nnode;
   long offset;
   sizs deltaPerimeter;

   list_to_propagate->curlev = lev;  /* Update queue curlev */

   deltaPerimeter = 0;

   nnode = tree->level[lev].nnodes + 1;

   while (!hfifo_empty(lev, list_to_propagate))
  {
    pixel  = gethqueue(list_to_propagate);   /* Get pixel from FIFO list */
    *(map-> f + (pixel - img->f)) = nnode;

    for(index = 0; index<connectivity; index++)
   {
     neighbour_pixel = pixel + direction[index];
     if (*neighbour_pixel == BORDER) continue;
     offset = neighbour_pixel - img->f;

     if (index<4){
       if (*neighbour_pixel<lev) deltaPerimeter++;
       else if (*neighbour_pixel>lev) deltaPerimeter--;
     }

     if (*(map->f+offset) != NOT_ANALYSED) continue; /* Already ANALYSED */

     node_at_level[*neighbour_pixel] = TRUE;

     *(map->f+offset) = IN_QUEUE;  /* Pixel in queue */
     putvqueue(neighbour_pixel, list_to_propagate);

     if (*neighbour_pixel>lev) /* The connected component has children */
    {
      l = *neighbour_pixel;

      do
     {
       /* Recursivity */
       flood(l, list_to_propagate, tree, map, img, connectivity, 
        direction, node_at_level, p_node);
       l = list_to_propagate->curlev;
     } 
      while (l>lev);  /* Check if there are any more 
          components with higher greylevel 
             than 'lev' */

    } /* if */

   } /* for (index = 0; index < 8; index ++) */

  } /* while (!fifo_pempty(lev, list_to_propagate))  */


   /* Allocate node */

   tree->level[lev].nnodes++;

   if ((nnode % EXPAND_SIZE_NODE) == 0) /* Do we have to allocate memory ? */
  {
    if (nnode == 0)
   tree->level[lev].node = (struct max_node *)
     emalloc((sizs) (sizeof(struct max_node) * EXPAND_SIZE_NODE));
    else
   tree->level[lev].node = (struct max_node *)
     realloc(tree->level[lev].node, 
       sizeof(struct max_node) * (nnode + EXPAND_SIZE_NODE));
    memset(tree->level[lev].node + nnode, 0x00, 
     sizeof(struct max_node) * EXPAND_SIZE_NODE );
  }

   p_node = tree->level[lev].node + nnode; /* point to actual node we study */

   p_node->RestitutionLevel = lev;
   p_node->deltaPerimeter   = deltaPerimeter;

   /* Search first non empty level of list */

   node_at_level[lev] = FALSE;
   do lev--; while ((lev>=0)&&(node_at_level[lev]==FALSE));
 
   if (lev >= 0)
  {
    p_node->father.level = lev;
    p_node->father.node  = tree->level[lev].nnodes + 1;
  }
   else
  {
    p_node->father.level = -1; /* We indicate that this node has no father */
    p_node->father.node  = -1;
  }
      
   list_to_propagate->curlev = lev;
      
 } /* void flood */


#undef NOT_ANALYSED
#undef IN_QUEUE
#undef BORDER


/******************************************************************************/
/*
 * -- compute_children --
 *
 * Purpose: This function computes the children relationships by knowing
 *  the father relationships (the latter ones are computed in create_max_tree).
 *  This function is called internally by create_max_tree. It should not be
 *  called externally.
 *
 * Arguments in/out:
 *  struct max_tree *tree
 *
 * Note: All information allocated in member "children" of "struct max_node"
 *  will be lost.
 *
 * Return values: None
 */
/***************************************************************************/

void compute_children (struct max_tree *tree)
{
    register datim lev, index;
    register struct max_node *p_node, *p_node_end, *p_dest;
    register struct point *p_father;
    
    /* Compute the number of children of each node */
    
    for(lev = 0; lev<GREY_LEVELS; lev++)  /* Scan all levels */
    {
        p_node     = tree->level[lev].node;
        p_node_end = p_node + tree->level[lev].nnodes;
        
        for(;p_node<p_node_end; p_node++)   /* Scan all nodes */
        {
            p_father = &(p_node->father);
            if (p_father->level != -1) /* Has p_node any father ? */
                tree->level[p_father->level].node[p_father->node].NChildren++;
        }
    }
    
    /* Allocate memory */
    
    for(lev = 0; lev<GREY_LEVELS; lev++)
    {
        p_node     = tree->level[lev].node;
        p_node_end = p_node + tree->level[lev].nnodes;
        
        for(; p_node<p_node_end; p_node++)
            if (p_node->NChildren != 0)
            {
                p_node->Child = (struct point *)
                emalloc((sizs) (sizeof(struct point) * p_node->NChildren));
                p_node->NChildren = 0; /* We'll use it as a counter */
            }
    }
    
    /* Allocate father info in children */
    
    for(lev = 0; lev<GREY_LEVELS; lev++)
    {
        p_node     = tree->level[lev].node;
        p_node_end = p_node + tree->level[lev].nnodes;
        
        index      = 0;
        
        for(; p_node<p_node_end; p_node++, index++)
        {
            p_father = &(p_node->father);
            if (p_father->level != -1)
            {
                p_dest   = tree->level[p_father->level].node + p_father->node;
                
                p_dest->Child[p_dest->NChildren].level = lev;
                p_dest->Child[p_dest->NChildren].node  = index;
                
                p_dest->NChildren++;
            }
        }
    }
}


/******************************************************************************/
/*
 *  undo_max_tree
 */
/*    Creates "out" (image representation) using "tree" (max tree
 *    representation). Image "out" should be previously allocated.
 *
 *  Once the max_tree has been filtered, the member
 *  RestitutionLevel of each node indicates to which gray level
 *  the pixels assigned to that node should be restituted. This 
 *  information is used to create "out".
 *
 */


void undo_max_tree (struct max_tree *tree, struct image *img)
{
  register datim lev;
  register sizs j;
  register struct max_node *p_node, *p_end;

  lev = -1;

  while (++lev < GREY_LEVELS) 
 {
      /* Restore the node with the corresponding grey level */

      p_node = tree->level[lev].node;  /* i.e. &(tree->level[lev].node[0])  */
      p_end  = p_node + tree->level[lev].nnodes;
      
      for(; p_node<p_end; p_node++) 
  { 
    for(j = 0; j<p_node->NPixels; j++)
            *(img->f + p_node->Pixel[j]) = p_node->RestitutionLevel;
  }
 }
}


/******************************************************************************/
/*
 *  free_max_tree
 */
/*    Free all the information allocated in max_tree. The
 *    same "struct max_tree" is not deallocated.
 *
 */


void free_max_tree(struct max_tree *tree)
{
  register datim lev;
  register struct max_node *p_node, *p_end;
   
  lev = -1;

  while (++lev<GREY_LEVELS) 
 {
      p_node = tree->level[lev].node;
      p_end  = p_node + tree->level[lev].nnodes;
      
      for(;p_node < p_end; p_node++) 
  {
    free((char *) p_node->Pixel);
    if (p_node->NChildren != 0) free((char *) p_node->Child);
  }
      
      if (tree->level[lev].nnodes!=0) free((char *) tree->level[lev].node);
      
 }

  free((char *) tree->level);
}


/**********************************************************************/
/* The next function is used internally by binary_connected_component */
/**********************************************************************/

long number_pixels_connected_component (struct max_node *p_node, 
          struct hqueue *list, 
          struct max_tree *tree)
{
     struct max_node *p_next_node;
     struct point *p_children, *p_children_end;
     long n_pixels;
     
     n_pixels = p_node->NPixels;

     putqueue((datim *) p_node, list);

     if (p_node->NChildren != 0)
     {
          p_children     = p_node->Child;
          p_children_end = p_children + p_node->NChildren;
     
          for(; p_children<p_children_end; p_children++)
          {
               p_next_node = tree->level[p_children->level].node + p_children->node;
               n_pixels += number_pixels_connected_component(p_next_node, list, tree);
          }
     }
     return n_pixels;
}

         
/******************************************************************************/
/*
 *  binary_connected_component
 */
/*    Computes the connected component associated to "node" of
 *    "tree". The result is returned in "pixel" and "npixels". Npixels is
 *    the number of pixels of the connected component, and pixel is a
 *    vector of "longs" containing the coordinates (the offsets) of the
 *    pixels belonging to this connected component. The vector pixel is
 *    allocated in this function.
 */
 

void binary_connected_component(struct max_tree *tree, struct max_node *p_node, 
        long **pixel, long *n_pixels)
{
     struct hqueue *list;

     list = new_fifo_queue();

     *n_pixels = number_pixels_connected_component(p_node, list, tree);

     *pixel = (long *) emalloc((sizs) (sizeof(long) * (*n_pixels)));

     *n_pixels = 0;
     
     while (!fifo_empty(list))
     {
          p_node = (struct max_node *) gethqueue(list);

          memcpy(*pixel + *n_pixels, p_node->Pixel, sizeof(long) * p_node->NPixels);
          *n_pixels += p_node->NPixels;
     }

     free_hqueue(list);
     free((char *) list);
}


/******************************************************************************/
/*
 *  Memory allocation and basic image related tools
 *
 */

void
runerr (    /*Run-time error handler.*/
        char *message,      /*Error message to be sent somewhere. (stderr)*/
        int exitcode        /*Program exit code.*/
)
{
    printf("%s \n", message);
    exit(exitcode);
}

char *
emalloc (   /*Memory allocation with error handling.
             WARNING: There is a problem when you want to allocate more memory
             than size can handle.  Malloc gets as argument an unsigned.  It poses
             a problem when sizeof(unsigned) < sizeof(char *)...*/
         sizs siz    /*Allocated chunk size.*/
)
{
    char *retval;
    if(siz <= 0) runerr("ERROR(emalloc): Wrong size: <= 0.",1);
    retval = malloc((unsigned int)siz);
    if(retval == NULL)
        runerr("ERROR(emalloc): memory allocation error.",1);
    return retval;
}

char *
ecalloc (   /*Memory allocation with error handling.
             WARNING: There is a problem when you want to allocate more memory
             than size can handle.  Malloc gets as argument an unsigned.  It poses
             a problem when sizeof(unsigned) < sizeof(char *)...*/
         sizs n,
         sizs s
         )
{
    char *retval;
    if(s <= 0 || n <= 0) runerr("ERROR(ecalloc): Wrong size: <= 0.",1);
    retval = calloc((unsigned int)n,(unsigned int)s);
    if(retval == NULL)
        runerr("ERROR(ecalloc): memory allocation error.",1);
    return retval;
}

void
imera ( /*Sets image data to a given value.*/
       struct image *a,    /*Input/output image.*/
       int val     /*Value.*/
)
{
    register datim *b,*e;
    e = a->f+a->x*a->y;
    for(b = a->f;b != e;b++) *b = val;
}

void
imfree (        /*Frees space provided by immalloc.*/
        struct image *o /*Image to be freed.*/
)
{
    if(o->x * o->y) free((char *)o->f);
    o->x = o->y = 0;
}

void
immalloc (  /*Image memory allocation.  Data is initialized to 0
             Structure itself must be allocated and initialized.*/
          struct image *out
          )
{
    if(out->grid != 's' && out->grid != 'h')
        runerr("ERROR(immalloc): invalid grid type",1);
    if(out->x > 0 && out->y > 0){
        out->version = VERSION;
        if (out->grid == 's') {
            out->upperodd = 0; /* Just to avoid Purify of complaining */
        }
        
        out->f = (datim *) ecalloc(out->x*out->y,(sizs)sizeof(datim));
        if(out->f == NULL) runerr("ERROR(immalloc): Memory allocation",1);
    }
}

void
timmalloc ( /*Image memory allocation for integer type, but
             with template templ.  It allocates the image according to the
             template templ.*/
           struct image *templ,    /*Template image.*/
           struct image *out       /*Image to be allocated.  Space for the
                                    structure itself must be allocated.*/
)
{
    out->x = templ->x;
    out->y = templ->y;
    out->grid = templ->grid;
    out->upperodd = templ->upperodd;
    immalloc(out);
    out->version = templ->version;
}

void
putbrd (    /*Puts borders on image.  Original image is modified.*/
        struct image *a,    /*Input/output image.*/
        int v           /*Value to be put on the border.*/
)
{
    datim *temp;
    register datim *s,*d;
    register sizs i,j;
    
    a->x += 2;
    a->y += 2;
    a->upperodd ^= 1;   /*Toggle upper line parity.*/
    
    temp = a->f;
    a->f = (datim *) emalloc((sizs) a->x*a->y*sizeof(datim));
    if(a->f == NULL) runerr("ERROR(putbrd): memory allocation error.",2);
    d = a->f;
    for(i = 0;i < a->x;i++) *d++ = v;
    s = temp;
    for(j = 1;j < a->y-1;j++){
        *d++ = v;
        for(i = 1;i != a->x-1;i++) *d++ = *s++;
        *d++ = v;
    }
    for(i = 0;i < a->x;i++) *d++ = v;
    free((char *)temp);
}

void
rmbrd ( /*Removes borders.*/
       struct image *a /*Input/output image.*/
)
{
    register datim *temp,*s,*d;
    register sizs i,j;
    
    a->x -= 2;
    a->y -= 2;
    a->upperodd ^= 1;   /*Toggle upper line parity.*/
    temp = a->f;
    a->f = (datim *) emalloc((sizs)a->x*a->y*sizeof(datim));
    if(a->f == NULL) runerr("ERROR(rmbrd): memory allocation error.",2);
    d = a->f;
    s = temp+a->x+3;
    for(j = 0;j != a->y;j++){
        for(i = 0;i != a->x;i++) *d++ = *s++;
        s += 2;
    }
    free((char *)temp);
}


/*****************************************************************************
 * Functions to handle fifo and hierarchial queues
 *
 */
/******************************************************************************/
/*
 *  expandhq
 */
/*   Expand the size of one of the queues involved in the hierarchical queue.
 *
 */

void
expandhq ( /*Expand a queue when there are too much points
            to fill it.*/
          struct hqueue *q, /*The queue.*/
          int level  /*The queue number to be expanded.  NOT the level of the
                      corresponding gradient image.*/
)
{
    register datim **r,**w,**b,**e,**nb;
    register sizs size;   /*Number of points needed for expansion.*/
    
    size = EXPAND_SIZE;
    /*Realloc an old queue.*/
    r = q->r[level];
    w = q->w[level];
    b = q->b[level];
    e = q->e[level];
    size += e - b;
    q->b[level] = (datim **) ecalloc(size + 1,(sizs) sizeof(datim *));
    q->e[level] = q->b[level] + size;
    nb = q->b[level];
    *nb = *r;
    nb++;
    r++;
    if(r >= e) r = b;
    while(r != w){
        *nb = *r;
        nb++;
        r++;
        if(r >= e) r = b;
    }
    q->r[level] = q->b[level];
    q->w[level] = nb;
    free((char *) b);
}

/******************************************************************************/
/*
 *  puthqueue
 */
/*  Put a pixel in the queue and mark it.
 */

void
puthqueue ( /*Put a pixel on a hierarchical queue.*/
           datim *val,
           struct hqueue *q
           )
{
    register datim level;
    
    if(*val == BIGVAL) return; /*Already on the queue...*/
    
    if(*val < q->curlev) level = q->curlev - q->minlev;
    else level = *val - q->minlev;
    /*This condition is there to perform homotopy
     modification on the flight.  Now, the queue number
     'level' is calculated.*/
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
    *val = BIGVAL; /*In order not to put it twice.*/
}

/******************************************************************************/
/*
 *  nputhqueue
 */
/*   Put a pixel in the queue without marking it.
 */
void
nputhqueue ( /*Put a pixel on a hierarchical queue.
              No flagging of the pixel.*/
            datim *val,
            struct hqueue *q
            )
{
    register datim level;
    
    if(*val == BIGVAL) return; /*Already on the queue...*/
    
    if(*val < q->curlev) level = q->curlev - q->minlev;
    else level = *val - q->minlev;
    /*This condition is there to perform homotopy
     modification on the flight.  Now, the queue number
     'level' is calculated.*/
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
}

/******************************************************************************/
/*
 *  nputhqueue_level
 */
/*  Put a pixel in the queue at a given level without marking it.
 *
 */

void
nputhqueue_level ( /*Put a pixel on a hierarchical queue.
                    No flagging of the pixel.*/
                  datim *val,
                  int lev,
                  struct hqueue *q
                  )
{
    register datim level;
    
    if(lev == BIGVAL) return; /*Already on the queue...*/
    
    if(lev < q->curlev) level = q->curlev - q->minlev;
    else level = lev - q->minlev;
    /*This condition is there to perform homotopy
     modification on the flight.  Now, the queue number
     'level' is calculated.*/
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
}


/******************************************************************************/
/*
 *  cputhqueue
 */
/*   Put a pixel in the queue at a given priority without marking it.
 *
 */

void
cputhqueue ( /*Put a pixel on a hierarchical queue.
              No flagging of the pixel.*/
            datim *val,
            int priority,
            struct hqueue *q
            )
{
    register datim level;
    
    
    if(priority < q->curlev) level = q->curlev - q->minlev;
    else level = priority - q->minlev;
    /*This condition is there to perform homotopy
     modification on the flight.  Now, the queue number
     'level' is calculated.*/
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
}

/******************************************************************************/
/*
 *  gethqueue
 */
/*   Extract the pixel of highest priority from the queue
 */

datim *
gethqueue ( /*Get a pixel on a hierarchical queue.*/
           struct hqueue *q
           )
{
    register datim maxlev,level;
    datim *oval;
    
    level = q->curlev - q->minlev;
    maxlev = q->maxlev - q->minlev;
    while(q->r[level] == q->w[level]){ /*Queue n empty.  Go to next queue.*/
        if(q->b[level] != NULL){ /*Free the queue.*/
            free((char *) q->b[level]);
            q->b[level] = q->r[level] = q->w[level] = NULL;
        }
        level++; /*Next level.*/
        if(level > maxlev){ /*Queue empty altogether.*/
            q->curlev = q->maxlev;
            return NULL;
        }
    } /*A non-empty queue is found.*/
    q->curlev = level + q->minlev; /*Curlev updated.*/
    
    oval = *(q->r[level]);
    q->r[level]++;
    if(q->r[level] >= q->e[level]){
        q->r[level] = q->b[level];
    }
    return oval;
}

/******************************************************************************/
/*
 *  kill_hqueue
 */
/*   Free memory allocated for a hierarchical queue.
 *
 */

void
kill_hqueue ( /*Free memory allocated for a hierarchical queue.*/
             struct hqueue *q
             )
{
    register datim level;
    
    for(level = 0; level <= q->maxlev + 1 - q->minlev;level++){
        if(q->b[level] != NULL){
            free((char *) q->b[level]);
            q->b[level] = NULL;
        }
    }
    free((char *)q->b);
    free((char *)q->e);
    free((char *)q->r);
    free((char *)q->w);
}

/******************************************************************************/
/*
 *  putvqueue
 */
/*    Put a pixel in the queue without marking it.
 *    The pixels with the lowest value have the highest priority.
 *    Arguments in:
 *    datim *val  Pixel to introduce in the queue and do not
 *    mark it
 */
void
putvqueue ( /*Put a pixel on a hierarchical queue without changing the value of
             the pixel. The pixels with the lowest value have the highest priority.*/
           datim *val,
           struct hqueue *q
           )
{
    register datim level;
    
    level = *val - q->minlev;
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
}

/******************************************************************************/
/*
 *  putpqueue
 */
/*    Put a pixel in the queue without marking it.
 *    The pixels with the highest value have the highest priority.
 *    Arguments in:
 *    datim *val  Pixel to introduce in the queue and do not
 *    mark it
 *
 */
void
putpqueue (  /*Put a pixel on a hierarchical queue without changing the value of
              the pixel. The pixels with the highest value have the highest priority.*/
           datim *val,
           struct hqueue *q
           )
{
    register datim level;
    
    level = q->maxlev - *val ;
    
    if(q->b[level] == NULL){ /*New subqueue.*/
        q->b[level]  = q->r[level] = q->w[level] = (datim **) ecalloc(
                                                                      (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[level] = q->b[level] + EXPAND_SIZE;
    }
    *(q->w[level]) = val;
    q->w[level]++;
    if(q->w[level] >= q->e[level]){
        q->w[level] = q->b[level];
    }
    if(q->w[level] == q->r[level]){
        expandhq(q,level);
    }
}


/******************************************************************************/
/*
 *  putqueue
 */
/*   Put a pixel in a first in fiirst out queue without marking it.
 *    This hierarchical queue acts as a simple fifo.
 */
void
putqueue (  /*Put a pixel on a queue FIRST IN/ FIRST OUT . Acts like a hierarchical queue with just one level. Same structure as for a hierarchical queue is used.*/
          datim *val,
          struct hqueue *q
          )
{
    
    if(q->b[0] == NULL){ /*New subqueue.*/
        q->b[0]  = q->r[0] = q->w[0] = (datim **) ecalloc(
                                                          (sizs)EXPAND_SIZE,(sizs)sizeof(datim *));
        q->e[0] = q->b[0] + EXPAND_SIZE;
    }
    *(q->w[0]) = val;
    q->w[0]++;
    if(q->w[0] >= q->e[0]){
        q->w[0] = q->b[0];
    }
    if(q->w[0] == q->r[0]){
        expandhq(q,0);
    }
}

/*****************************************************************************/
/*
 *  new_fifo_hqueue
 */
/*   Initialization of hierarchical queue
 *
 */

struct hqueue *new_fifo_hqueue(int max_val, int min_val)
{
    register struct hqueue *queue;
    register int nlevels;
    
    //queue = (struct hqueue *) emalloc( (sizs) sizeof(struct hqueue));
    queue = (struct hqueue *) emalloc( sizeof(struct hqueue));
    
    nlevels = max_val - min_val + 2;
    
    //  queue->b = (datim ***) ecalloc((sizs) nlevels, (sizs) sizeof(datim **));
    //  queue->e = (datim ***) ecalloc((sizs) nlevels, (sizs) sizeof(datim **));
    //  queue->r = (datim ***) ecalloc((sizs) nlevels, (sizs) sizeof(datim **));
    //  queue->w = (datim ***) ecalloc((sizs) nlevels, (sizs) sizeof(datim **));
    queue->b = (datim ***) ecalloc(nlevels,  sizeof(datim **));
    queue->e = (datim ***) ecalloc(nlevels,  sizeof(datim **));
    queue->r = (datim ***) ecalloc(nlevels,  sizeof(datim **));
    queue->w = (datim ***) ecalloc(nlevels,  sizeof(datim **));
    queue->maxlev = max_val;
    queue->minlev = queue->curlev = min_val;
    
    return queue;
}


/*****************************************************************************/
/*
 *  nitems_fifo_hqueue
 */
/*    Returnt the number of elements of given priority
 *    in a hierachical queue
 *
 */

long nitems_fifo_hqueue(struct hqueue *queue, int level)
{
    register long num;
    
    num = queue->w[level] - queue->r[level];
    if (num < 0) num = queue->e[level] - queue->b[level] + num;
    return num;
}




