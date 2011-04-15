#ifndef _NEUZEICHNEN_H
#define _NEUZEICHNEN_H

#include "basicdef.h"
#include "menu.h"

#define colorscalenum 3
extern Colormap     cmap;


/* Definition of differnt types of graphic routines that need 
   different types of data for definerequests.*/
#define curve_2d     0
#define contour_u_2d 1
#define graph_2d     2
#define contour_3d   3
#define project_3d   4

extern void SetComment(char *newcomment);

extern void definerequests(grstr *graph,int type,nzdata *zd);

extern void neu_zeichnen_1(Display *myd,mywindow *mw,
			   int PS,char *psname, double PSwidth, 
			   double PSheight, int PSletter, int PSlandscape,
			   grstr);
extern void neu_zeichnen_grfft(Display *myd,mywindow *mw,
			       int PS,char *psname, double PSwidth, 
			       double PSheight, int PSletter, int PSlandscape,
			       grstr);
extern void neu_zeichnen_Structfact(Display *myd,mywindow *mw,
				    int PS,char *psname, double PSwidth, 
				    double PSheight, int PSletter, int PSlandscape,
				    grstr);
extern void neu_zeichnen2(Display *myd,mywindow *mw,
			  int PS,char *psname, double PSwidth, 
			  double PSheight, int PSletter, int PSlandscape,
			  grstr);
extern void neu_zeichnen_contour3d(Display *myd,mywindow *mw,
				   int PS,char *psname, double PSwidth, 
				   double PSheight, int PSletter, int PSlandscape,
				   grstr);
extern void neu_zeichnen_data(Display *myd,mywindow *mw,
			      int PS,char *psname, double PSwidth, 
			      double PSheight, int PSletter, int PSlandscape,
			      grstr);

extern void neu_zeichnen_time(Display *myd,mywindow *mw,
			      int PS,char *psname, double PSwidth, 
			      double PSheight, int PSletter, int PSlandscape,
			      grstr);
extern void neu_zeichnen_2dgraph(Display *myd,mywindow *mw,
				 int PS,char *psname, double PSwidth, 
				 double PSheight, int PSletter, int PSlandscape,
				 grstr graph);

int Contour_print(char *psname, double PSwidth, double PSheight, 
                   int PSletter, int PSlandscape,
		   int nocuts,  double *cut, double *gp, int nox, int noy,
		   int comments, const char *comment, int magnify, 
		   int center,
		   int NewXsize, int NewYsize,double xofs, double yofs,
		   int adjustcuts,
		   int adjustdensity, double density_min, 
		   double density_max, int density);
#endif
