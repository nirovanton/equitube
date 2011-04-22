#ifndef GRAPH_H
#define GRAPH_H

#include "basicdef.h"

#define N_R         0
#define N_RxR       1
#define N_RxRxR     2
#define N_RxRxRxR   12
#define NxN_R       3
#define NxN_RxR     4
#define NxN_RxRxRt  5
#define cont3d      6
#define N_RxRp      7
#define N_Rp        8
#define NxNxN_R     9
#define NxNxN_Rp    10
#define NxN_Rp      11

typedef int *pint; 
typedef pint dim1[1];
typedef pint dim2[2];
typedef pint dim3[3]; 

typedef struct _grdat {
  int color,linetype,shape,fill;
  double size;
} grdat;

typedef struct _grstr {  
  int anzN_R,maxanzN_R,          /*2d-curve*/
      anzN_Rp,maxanzN_Rp,          /*2d-curve*/
      anzN_RxR,maxanzN_RxR,      /*2d-curve*/
      anzN_RxRp,maxanzN_RxRp,      /*2d-curve*/
      anzN_RxRxR,maxanzN_RxRxR,  /*3-d curve, not implemented*/
      anzN_RxRxRxR,maxanzN_RxRxRxR,  /*2-d Tubes for John Harris*/
      anzNxN_R,maxanzNxN_R,      /*3-d graph, 2d-contourplot*/
      anzNxN_Rp,maxanzNxN_Rp,      /*3-d graph, 2d-contourplot*/
      anzNxN_RxR,maxanzNxN_RxR,  /*2-d vectorplot*/
      anzNxN_RxRp,maxanzNxN_RxRp,  /*2-d vectorplot*/
      anzNxN_RxRxRt,maxanzNxN_RxRxRt,  /*2-d symmetric tensor*/
      anzNxNxN_R,maxanzNxNxN_R,  /*3-d plots of density/iterface*/
      anzNxNxN_Rp,maxanzNxNxN_Rp,  /*3-d plots of density/iterface*/
      anzNxNxN_RxRxR,maxanzNxNxN_RxRxR,/*not implemented*/
      anzcont3d,maxanzcont3d;
  char **N_Rname,**N_Rnamep,**N_RxRname,**N_RxRpname,**N_RxRxRname,**N_RxRxRxRname,
       **NxN_Rname,**NxN_Rpname,**NxN_RxRname,**NxN_RxRpname,**NxN_RxRxRtname,
       **NxNxN_Rname,**NxNxN_Rnamep,**NxNxN_RxRname,
       **cont3dname;
  double **pN_R,***pN_Rp,**pN_RxR,***pN_RxRp,**pN_RxRxR,**pN_RxRxRxR,
    **pNxN_R,***pNxN_Rp,**pNxN_RxR,***pNxN_RxRp,**pNxN_RxRxRt,
    **pNxN_RxRxR,**pNxNxN_R,***pNxNxN_Rp,**pNxNxN_RxRxR;
  grdat *dN_R,*dN_Rp,*dN_RxR,*dN_RxRp,*dN_RxRxRxR; /* only some types require these data */
  Contour ***pcont3d;
  dim1    *dimN_R,*dimN_Rp,*dimN_RxR,*dimN_RxRp,*dimN_RxRxR,*dimN_RxRxRxR;
  dim2    *dimNxN_R,*dimNxN_Rp,*dimNxN_RxR,*dimNxN_RxRp,*dimNxN_RxRxRt;
  dim3    *dimNxNxN_R,*dimNxNxN_Rp,*dimNxNxN_RxRxR,*dimcont3d;
  int **reqN_R,**reqN_Rp,**reqN_RxR,**reqN_RxRp,**reqN_RxRxR,**reqN_RxRxRxR,**reqNxN_R,**reqNxN_Rp,
      **reqNxN_RxR,**reqNxN_RxRp,**reqNxN_RxRxRt,**reqNxN_RxRxR,
      **reqNxNxN_R,**reqNxNxN_Rp,**reqNxNxN_RxRxR,**reqcont3d;
  /* This was a nice idea that worked beautifully with c-programs.
     Since ZPL does not support functional parameters we are required
     to move to a less general and satisfactory solution. 
  gr1func *getN_R,*getN_RxR,*getN_RxRxR;
  gr2func *getNxN_R,*getNxN_RxR,*getNxN_RxRp;
  gr3func *getNxNxN_R,*getNxNxN_RxRxR;
  */

} grstr; 

extern grstr *graph;
extern int graphactive;

extern void SetDefaultColor(int c);
extern void SetDefaultLineType(int s);
extern void SetDefaultShape(int s);
extern void SetDefaultSize(double s);
extern void SetDefaultFill(int f);


extern void NewGraph();
extern void SetActiveGraph(int);

/*extern void DefineGraph(grstr *graph,char *Name,int Type,...);*/
extern void DefineGraphN_R(char *Name,double *gd,int *dim,int *req);
extern void DefineGraphN_Rp(char *Name,double **gd,int *dim,int *req);
extern void DefineGraphN_RxR(char *Name,double *gd,int *dim,int *req);
extern void DefineGraphN_RxRxRxR(char *Name,double *gd,int *dim,int *req);
extern void DefineGraphN_RxRp(char *Name,double **gd,int *dim,int *req);
extern void DefineGraphNxN_R(char *Name,double *gd,
                      int *dim1,int *dim2,int *req);
extern void DefineGraphNxN_Rp(char *Name,double **gd,
                      int *dim1,int *dim2,int *req);
extern void DefineGraphNxN_RxR(char *Name,double *gd,
                      int *dim1,int *dim2,int *req);
extern void DefineGraphNxN_RxRxRt(char *Name,double *gd,
                      int *dim1,int *dim2,int *req);
extern void DefineGraphContour3d(char *Name,Contour **gd,
                      int *dim1,int *dim2,int *dim3,int *req);
extern void getgraph3d(grstr graph,int type,double **gp,int *no,int get, 
		int rough, int *copy);
#endif
