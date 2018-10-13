/* ==========================================================
 * UNIVERSITAT POLITECNICA DE CATALUNYA, Spain
 * Image processing group: https://imatge.upc.edu/web/
 * ==========================================================
 *
 * max_tree.h
 * Include for the MATLAB mex Interface maxtree_creation.c
 *
 * Compute the maxtree representation on a matlab int16 image 
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
 *		maxtree = maxtree_creation(image, conectivity)
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

/******************        INCLUDE FILES       *******************************/
#   include <stdlib.h>                       /* SURE ALL Def. ARE DEFINED !  */
#   include <math.h>                         /********************************/
#   include <string.h>    

/* ================================================== */

#   if !defined(FALSE) || ((FALSE)!= 0)      /* TO AVOID MULTIPLE DEFINE     */
#      define  FALSE           0             /* AND BE SURE FALSE = 0        */
#   endif                                    /********************************/
#   if !defined(TRUE)  || ((TRUE) != 1)      /* TO AVOID MULTIPLE DEFINE     */
#      define  TRUE            1             /* AND BE SURE TRUE = 1         */
#   endif                                    /********************************/

#   ifndef     NULL
#      define  NULL            0
#   endif


/* ================================================================== */
/*                  Image structure and access                        */
typedef long  sizs;         /*Data size type.  MUST BE SIGNED!*/
typedef short int datim;    /*Image data type.*/

struct image{
    long version;  /*!< Version number. */
    sizs x;        /*!< Image width in pixels. */
    sizs y;        /*!< Image height in pixels. */
    char upperodd; /*!< Flag to tell if the top line is considered as
                    even or odd.  This is used for hex grids.*/
    char grid;     /*!< Grid type: s = square, h = hexagonal.*/
    datim *f;
};

#define PIX(m, i, j) ((m)->f[(j) * (m)->x + (i)])
#define PIX_COOR_Y(m,p) ((p) / (m)->x)
#define PIX_COOR_X(m,p) ((p) % (m)->x)
#define IMAGE_X(m) ((m)->x)
#define IMAGE_Y(m) ((m)->y)

/* ================================================================== */
#define VERSION 2 /*Lib version.*/

#define BIGVAL 32767                       /*Max value of short image data.*/
#define SMLVAL -32767                      /*A small value for image data.*/

/* Number of image gray-levels. Any positive number ( < 32001) should work */
#define GREY_LEVELS         32001

/* ================================================================== */
/*      Declaration related to fifo and hierarchical queues           */

struct hqueue{  /*Structure for hierarchical queues.*/
    datim ***b,***e; /*Begin and end of all the queues.*/
    datim ***r,***w; /*Read and write pointers of the queues.*/
    datim curlev,minlev,maxlev;
    /*Different levels in the image: current, minimum, maximum.*/
};

struct hqueue *new_fifo_hqueue(int max_val, int min_val);
long nitems_fifo_hqueue(struct hqueue *queue, int level);

/* Init a non hierarchical FIFO queue */
#define new_fifo_queue()                    new_fifo_hqueue(0, 0)

/* Check if FIFO is empty */
#define fifo_empty(q)                       ((q)->r[0]   == (q)->w[0])
#define hfifo_empty(k, q)                   ((q)->r[(k)] == (q)->w[(k)])

/* Flush FIFO queues */
#define flush_fifo_queue(q)                 ((q)->r[0]   = (q)->w[0] = (q)->b[0])
#define flush_hfifo_queue(k, q)             ((q)->r[(k)] = (q)->w[(k)])

/* Compute the number of items of the FIFO queue */
#define nitems_fifo_queue(q)                nitems_fifo_hqueue((q), 0)
#define nitems_hfifo_queue(k, q)            nitems_fifo_hqueue((q), (k))

/* Expand_size_node is used to create the max tree. A larger number
   will speed up the max tree creation, but needs more memory. A
   small number will decrease a lot the creation */
#define EXPAND_SIZE_NODE    100
/* Size by which a queue will be expanded if there is not enough memory.  
   Larger will be faster, but uses up more memory.*/
#define EXPAND_SIZE         3000

#define free_hqueue(q)      kill_hqueue(q)

/* Point is used to specify a tree node in the matrix structure where it is allocated */
struct point {
   datim level;
   datim node;
};

/* The node structure */
struct max_node {
   long  nodeID;
   struct point father;

   datim NChildren;
   struct point *Child;

   long NPixels;  /* Number of pixels associated to the node */
   long *Pixel;   /* A pixel is allocated as an offset from base address */
  
   datim RestitutionLevel;
   sizs deltaPerimeter;
};

/* An array of max tree nodes */
struct array_max_nodes {
   datim nnodes;                    /* Number of nodes */
   struct max_node *node;           /* Array of nodes */
};

/* And finally the max tree structure */
struct max_tree {
   sizs image_sizx, image_sizy;     /* Extra data. Size of the image. */
   struct array_max_nodes *level;   /* Pointer to the tree */
};

/* Functions declarations */
void runerr(char *message, int exitcode);
char *emalloc(sizs siz);
char *ecalloc(sizs n, sizs s);
void immalloc(struct image *out);
void timmalloc(struct image *templ, struct image *out);
void imfree(struct image *o);
void imera(struct image *a, int val);
void ran(struct image *a, int amp);
void imcopy(struct image *a, struct image **o);
void imcp(struct image *a, struct image *o);
void chkim(struct image *a, struct image *b, char *from, int grid);
void notsame(struct image *a, struct image *b, char *from);
void putbrd(struct image *a, int v);
void rmbrd(struct image *a);

void expandhq(struct hqueue *q, int level);
void puthqueue(datim *val, struct hqueue *q);
void nputhqueue(datim *val, struct hqueue *q);
void nputhqueue_level(datim *val, int lev, struct hqueue *q);
void cputhqueue(datim *val, int priority, struct hqueue *q);
datim *gethqueue(struct hqueue *q);
void kill_hqueue(struct hqueue *q);
void putvqueue(datim *val, struct hqueue *q);
void putpqueue(datim *val, struct hqueue *q);
void putqueue(datim *val, struct hqueue *q);
struct hqueue *new_fifo_hqueue(int max_val, int min_val);
long nitems_fifo_hqueue(struct hqueue *queue, int level);

void compute_children(struct max_tree *tree);
void flood(datim lev, struct hqueue *list_to_propagate, struct max_tree *tree, struct image *map, struct image *img, int connectivity, int *direction, datim *node_at_level, struct max_node *p_node);
void create_max_tree(struct max_tree *tree, struct image *img, int connectivity);
void undo_max_tree(struct max_tree *tree, struct image *img);
void free_max_tree(struct max_tree *tree);
long number_pixels_connected_component(struct max_node *p_node, struct hqueue *list, struct max_tree *tree);
void binary_connected_component(struct max_tree *tree, struct max_node *p_node, long **pixel, long *n_pixels);



