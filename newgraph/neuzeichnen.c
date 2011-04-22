/***************************************************************************\
graph library
Copyright (C) 1992-2010  Alexander Wagner

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

e-mail: Alexander.Wagner@ndsu.nodak.edu
\***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mypow.h"
#include "mydraw.h"
#include "mygraph.h"
#include "cobject3d.h"
#include "koord2d.h"
#include "draw2dgraph.h"
#include "draw3dgraph.h"
#include "draw3dcontour.h"
#include "ufield.h"
#include "tfield.h"
#include "contour.h"
#include "neuzeichnen.h"

static char *UserComment=NULL;

void SetComment(char *newcomment){
  UserComment=newcomment;
}

void definerequests(grstr *graph,int type,nzdata *zd){
  int i,j;

  switch (type){
  case curve_2d:
      for (i=0;i<graph->anzN_R;i++)
       if (graph->reqN_R[i]!=NULL)
	 *graph->reqN_R[i]|=zd->graph2d.draw1[i];
      for (j=0;j<graph->anzN_Rp;i++,j++)
       if (graph->reqN_Rp[j]!=NULL)
	 *graph->reqN_Rp[j]|=zd->graph2d.draw1[i];
      for (i=0;i<graph->anzN_RxR;i++)
       if (graph->reqN_RxR[i]!=NULL)
	 *graph->reqN_RxR[i]|=zd->graph2d.draw[i];
      for (j=0;j<graph->anzN_RxRp;i++,j++)
       if (graph->reqN_RxRp[j]!=NULL)
	 *graph->reqN_RxRp[j]|=zd->graph2d.draw[i];
      for (i=0;i<graph->anzN_RxRxRxR;i++)
       if (graph->reqN_RxRxRxR[i]!=NULL)
	 *graph->reqN_RxRxRxR[i]|=zd->graph2d.draw2[i];
    break;
  case graph_2d:
    for (i=0;i<graph->anzNxN_R;i++)
       if (graph->reqNxN_R[i]!=NULL)
	 *graph->reqNxN_R[i]|=zd->nz1.getgraph[i];
    for (i=0;i<graph->anzNxN_Rp;i++)
       if (graph->reqNxN_Rp[i]!=NULL)
	 *graph->reqNxN_Rp[i]|=zd->nz1.getgraph[i];
    zd->nz1.newdata=1;
    break;
  case contour_u_2d:
    for (i=0;i<graph->anzNxN_R;i++)
       if (graph->reqNxN_R[i])
	 *graph->reqNxN_R[i]|=zd->nz2.draw[i];
    for (i=0;i<graph->anzNxN_Rp;i++)
       if (graph->reqNxN_Rp[i])
	 *graph->reqNxN_Rp[i]|=zd->nz2.draw[i+graph->anzNxN_R];
    for (i=0;i<graph->anzNxN_RxR;i++)
       if (graph->reqNxN_RxR[i])
	 *graph->reqNxN_RxR[i]|=zd->nz2.draw_ufield[i];
    for (i=0;i<graph->anzNxN_RxRxRt;i++)
       if (graph->reqNxN_RxRxRt[i])
	 *graph->reqNxN_RxRxRt[i]|=zd->nz2.draw_tfield[i];
    break;
  case project_3d:
    for (i=0;i<graph->anzNxNxN_R;i++)
       if (graph->reqNxNxN_R[i])
	 *graph->reqNxNxN_R[i]|=
	   zd->p3d.drawX[i]||zd->p3d.drawY[i]||zd->p3d.drawZ[i];

    break;

  case contour_3d:
    for (i=0;i<graph->anzcont3d;i++)
       if (graph->reqcont3d[i])
	 *graph->reqcont3d[i]|=zd->nzc.getgraph[i];
    zd->nzc.newdata=1;
    break;
  default:
    printf("Error in neuzeichnen.c:definerequest: type %i unknown.\n",
	   type);
  }
}

void scalegraph(double **gp,int *getgraph,int anz,long no,
		double *zmin, double *zmax)
{
  int i,init=0;
  long j;
  double dz,*lgp;

  for (i=0;(getgraph[i]==0)&&(i<anz);i++);
  if (i==anz) return;
  for (; i<anz; i++) 
    {
      if (getgraph[i]!=0)
	{
	  lgp=gp[i]-1;
	  for (j=0; j<no; j++)
	    {
	      if (init==0){
		if (isnormal(*(lgp+1))){
		  *zmin = *zmax = *(lgp+1);
		  init=1;
		}
	      }
	      else {
		if (*(++lgp) < *zmin) 
		  if (isnormal(*lgp)) *zmin = *lgp;
		if (*lgp > *zmax) 
		  if (isnormal(*lgp)) *zmax = *lgp;
	      }
	    }
	}
    }
  dz=(*zmax-*zmin)/50;
  if (dz==0) dz=1;
  *zmax+=dz;
  *zmin-=dz;
}

void scaleScalar(double *gp,long no,double *zmin, double *zmax)
{
  long j;
  int init=0;
  double dz,*lgp;

  lgp=gp-1;
  for (j=0; j<no; j++)
    {
      if (init==0){
	if (isnormal(*(++lgp))){
	  *zmin = *zmax = *lgp;
	  init=1;
	}
      } else {
	if (*(++lgp) < *zmin) 
	  if (isnormal(*lgp)) *zmin = *lgp;
	if (*lgp > *zmax) 
	  if (isnormal(*lgp)) *zmax = *lgp;
      }
      }
  dz=(*zmax-*zmin);
  if (dz==0){
    dz=1;
    *zmax+=dz;
    *zmin-=dz;
  }
}

void scalegraph1d(double **gp,int *getgraph,int anz,dim1 *dim,
		  double *xmin, double *xmax,double *ymin, double *ymax)
{
  int i;
  long j;
  double *lgp;

  /* This algorithm assumes that xmin,xmax,ymin and ymax have been
     set to appropriate values already.*/
    
  for (i=0;(getgraph[i]==0)&&(i<anz);i++);
  if (i==anz) return;
  if ((*xmin==0) && (*xmax==0)&&(*ymin==0)&&(*ymax==0)){
    *xmin = *xmax = 0;
    *ymin = *ymax = gp[i][1];
  }
  for (; i<anz; i++) {
    if (getgraph[i]!=0){
      lgp=gp[i]-1;
      if (*xmin>0) *xmin=0;
      if (*xmax<*dim[i][0]) *xmax=*dim[i][0];
      for (j=0; j<*dim[i][0]; j++){
	if (*(++lgp) < *ymin) *ymin = *lgp;
	if (   *lgp  > *ymax) *ymax = *lgp;
      }
    }
  }
}

void scalegraph2d(double **gp,int *draw,int anz,
		 dim1 *dim, double *xmin,
		 double *xmax,double *ymin,
		 double *ymax)
{
  int i,count;
  long j;
  double *lgp;

  // is there a drawable graph?
  for (i=0,count=0;i<anz;i++){
    if ((draw[i]==0)||(*dim[i][0]==0)) count++;
    else i=anz;
  }
  if (count==anz) return;
  i=count;

  *xmin = *xmax = gp[i][0];
  *ymin = *ymax = gp[i][1];
  for (; i<anz; i++) 
    {
      if ((draw[i]!=0)&&(*dim[i][0]!=0))
	{
	  lgp=gp[i]-1;
	  for (j=0; j<*dim[i][0]; j++)
	    {
	      if (*(++lgp) < *xmin) 
		*xmin = *lgp;
	      if (*lgp > *xmax) 
		*xmax = *lgp;
	      if (*(++lgp) < *ymin) 
		*ymin = *lgp;
	      if (*lgp > *ymax) 
		*ymax = *lgp;
	    }
	}
    }
}

void scalegraph_tubes(double **gp,int *draw,int anz,
		 dim1 *dim, double *xmin,
		 double *xmax,double *ymin,
		 double *ymax)
{
  int i,count;
  long j;
  double *lgp;

  // is there a drawable graph?
  for (i=0,count=0;i<anz;i++){
    if ((draw[i]==0)||(*dim[i][0]==0)) count++;
    else i=anz;
  }
  if (count==anz) return;
  i=count;

  *xmin = *xmax = gp[i][0];
  *ymin = *ymax = gp[i][1];
  for (; i<anz; i++) 
    {
      if ((draw[i]!=0)&&(*dim[i][0]!=0))
	{
	  lgp=gp[i]-1;
	  for (j=0; j<*dim[i][0]; j++)
	    {
	      if (*(++lgp) < *xmin) 
		*xmin = *lgp;
	      if (*lgp > *xmax) 
		*xmax = *lgp;
	      if (*(++lgp) < *ymin) 
		*ymin = *lgp;
	      if (*lgp > *ymax) 
		*ymax = *lgp;
	      lgp+=2; /* angle and length */
	    }
	}
    }
}

void AddBorder(double *xmin,double *xmax,double *ymin,double *ymax){
  double dz;

  dz=(*xmax-*xmin)/50;
  if (dz==0) dz=1;
  *xmax+=dz;
  *xmin-=dz;
  dz=(*ymax-*ymin)/50;
  if (dz==0) dz=1;
  *ymax+=dz;
  *ymin-=dz;
}

double vec_scale(int no,double *vp)
{
  register int i;
  register double scale,tmp;
  
  scale=0;
  for (i=0;i<2*no;i+=2){
    tmp=pow2(vp[i])+pow2(vp[i+1]);
    if (scale<tmp) scale=tmp;
  }
  return sqrt(scale);
}

double te_scale(int no,double *vp)
{
  register int i;
  register double scale,tmp;
  
  scale=0;
  for (i=0;i<3*no;i+=3){
    tmp=fabs(vp[i])+fabs(vp[i+1])+1.42*fabs(vp[i+2]);
    if (scale<tmp) scale=tmp;
  }
  return scale*1.1;
}

double vec_scale_mask(int no, double *vp, double *zp, double cut)
{
  register int i;
  register double scale,tmp;

  scale=0;
  for (i=0;i<no;i++)
      {
	if (zp[i]<cut) {vp[2*i]=vp[2*i+1]=0;}
        else{
	  tmp=sqrt(pow2(vp[2*i])+pow2(vp[2*i+1]));
	  if (scale<tmp)
	    scale=tmp;
	}
      }
  return scale;
}

void neu_zeichnen_1(Display *myd,mywindow *mw,int PS,char *psname, double PSwidth, double PSheight, int PSletter, int PSlandscape,
		   grstr graph)
{
  static int i,xsize,ysize,no[2],xmin,xmax,ymin,ymax;
  static int initialized=0;
  static double **gp;
  int graphexist=0;
  
  if (initialized==0){
    gp=(double **) malloc((graph.anzNxN_R+graph.anzNxN_Rp)*sizeof(double *));
    for (i=0;i<graph.anzNxN_R+graph.anzNxN_Rp;i++) gp[i]=NULL;
    initialized=1;
  }
  /*  if (mw->menu->zd->nz1.newdata)
    {*/
      for (i=0;i<graph.anzNxN_R;i++)
	if (mw->menu->zd->nz1.getgraph[i]){
	  graphexist=1;
	  getgraph3d(graph,NxN_R,&(gp[i]),
		 &no[0],i,mw->menu->zd->nz1.rough,
		     &mw->menu->zd->nz1.copy[0]);
	}
      for (i=0;i<graph.anzNxN_Rp;i++)
	if (mw->menu->zd->nz1.getgraph[graph.anzNxN_R+i]){
	  graphexist=1;
	  getgraph3d(graph,NxN_R,&(gp[graph.anzNxN_R+i]),
		 &no[0],graph.anzNxN_R+i,mw->menu->zd->nz1.rough,
		     &mw->menu->zd->nz1.copy[0]);
	}
      /*if (graphexist==0) return; 
      necessary because of missing dimensions for initialization of mesh*/
      if (mw->menu->zd->nz1.Autoscale||mw->menu->zd->nz1.scaleonce)
	scalegraph(gp,mw->menu->zd->nz1.getgraph,graph.anzNxN_R,
		   no[0]*no[1],
		   &mw->menu->zd->nz1.zmin,&mw->menu->zd->nz1.zmax);
      /* should add scaling for NxN_Rp data type */
      mw->menu->zd->nz1.scaleonce=0;
      /*    }*/
  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }

  xmin=0; xmax=no[0]*mw->menu->zd->nz1.rough; 
  ymin=0; ymax=no[1]*mw->menu->zd->nz1.rough;
  draw3dgraph(xsize,ysize,graph,gp,mw->menu->zd->nz1.getgraph,
	      no[0],no[1],
	      &mw->menu->zd->nz1.NEWDIMENSIONS,mw->menu->zd->nz1.newdata,
	      /*mw->menu->zd->nz1.comments*/0,
	      /*get_time()*/ 0,
	      mw->menu->zd->nz1.persp_factor,mw->menu->zd->nz1.scalefakt, 
	      mw->menu->zd->nz1.xangle,mw->menu->zd->nz1.yangle,0,
	      xmax,xmin,ymax,ymin,
	      mw->menu->zd->nz1.zmin,mw->menu->zd->nz1.zmax,
	      mw->menu->zd->nz1.zfact,
	      &(mw->menu->zd->nz1.me),&(mw->menu->zd->nz1.meinit));
  mw->menu->zd->nz1.newdata=0;
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}


void neu_zeichnen_2dgraph(Display *myd,mywindow *mw,int PS,char *psname, double PSwidth, double PSheight, int PSletter, int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,i,j,done=0;
  dim1 *dimN_R,*dimN_RxR,*dimN_RxRxRxR;
  static  double **gp=NULL,**gp1=NULL,**gp2=NULL;
  grdat *mygrdat,*mygrdat1;

  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }

  mygrdat=(grdat *)malloc((graph.anzN_RxR+graph.anzN_RxRp)*sizeof(grdat));
  memcpy(mygrdat,graph.dN_RxR,graph.anzN_RxR*sizeof(grdat));
  memcpy(&mygrdat[graph.anzN_RxR],graph.dN_RxRp,graph.anzN_RxRp*sizeof(grdat));

  gp= (double **) malloc((graph.anzN_RxR+graph.anzN_RxRp)*sizeof(double *));
  if (gp==NULL) {
    fprintf(stderr,"neuzeichen.c:2dgraph:can not allocate gp:%i*sizeof(double*)\n",
	    graph.anzN_RxR+graph.anzN_RxR);
    exit(1);
  }
  mygrdat1=(grdat *)malloc((graph.anzN_R+graph.anzN_Rp)*sizeof(grdat));
  memcpy(mygrdat1,graph.dN_R,graph.anzN_R*sizeof(grdat));
  memcpy(&mygrdat1[graph.anzN_R],graph.dN_Rp,graph.anzN_Rp*sizeof(grdat));

  gp1= (double **) malloc((graph.anzN_R+graph.anzN_Rp)*sizeof(double *));
  if (gp1==NULL) {
    fprintf(stderr,"neuzeichen.c:2dgraph:can not allocate gp1:%i*sizeof(double*)\n",
	    graph.anzN_R+graph.anzN_Rp);
    exit(1);
  }

  gp2= (double **) malloc((graph.anzN_RxRxRxR)*sizeof(double *));
  if (gp1==NULL) {
    fprintf(stderr,"neuzeichen.c:2dgraph:can not allocate gp2:%i*sizeof(double*)\n",
	    graph.anzN_RxRxRxR);
    exit(1);
  }

  for (i=0; i< graph.anzN_RxR+graph.anzN_RxRp; i++){
    gp[i]=NULL; /* Why this? */
    if (mw->menu->zd->graph2d.draw[i]){
      done=1;
    }
  }
  for (i=0; i< graph.anzN_R+graph.anzN_Rp; i++){
    gp1[i]=NULL; /* Why this? */
    if (mw->menu->zd->graph2d.draw1[i]){
      done=1;
    }
  }
  for (i=0; i< graph.anzN_RxRxRxR; i++){
    gp2[i]=NULL; /* Why this? */
    if (mw->menu->zd->graph2d.draw2[i]){
      done=1;
    }
  }
  if (done==0) {
    free(mygrdat);
    free(mygrdat1);
    if (gp!=NULL) free(gp);
    if (gp1!=NULL) free(gp1);
    if (gp2!=NULL) free(gp2);
    myclear();
    myshow();
    return;
  }


  for (i=0; i< graph.anzN_R; i++){
    if (mw->menu->zd->graph2d.draw1[i])
      getgraph3d(graph,N_R,&gp1[i], graph.dimN_R[i][0],i,
		 mw->menu->zd->graph2d.rough,
		 &mw->menu->zd->graph2d.copy);
  }
  for (j=0; j< graph.anzN_Rp; i++,j++){
    if (mw->menu->zd->graph2d.draw1[i])
      getgraph3d(graph,N_Rp,&gp1[i], graph.dimN_Rp[j][0],j,
		 mw->menu->zd->graph2d.rough,
		 &mw->menu->zd->graph2d.copy);
  }
  for (i=0; i< graph.anzN_RxR; i++){
    if (mw->menu->zd->graph2d.draw[i])
      getgraph3d(graph,N_RxR,&gp[i], graph.dimN_RxR[i][0],i,
		 mw->menu->zd->graph2d.rough,
		 &mw->menu->zd->graph2d.copy);
  }
  for (j=0; j< graph.anzN_RxRp; i++,j++){
    if (mw->menu->zd->graph2d.draw[i])
      getgraph3d(graph,N_RxRp,&gp[i], graph.dimN_RxRp[j][0],j,
		 mw->menu->zd->graph2d.rough,
		 &mw->menu->zd->graph2d.copy);
  }
  for (i=0; i< graph.anzN_RxRxRxR; i++){
    if (mw->menu->zd->graph2d.draw2[i])
      getgraph3d(graph,N_RxRxRxR,&gp2[i], graph.dimN_RxRxRxR[i][0],i,
		 mw->menu->zd->graph2d.rough,
		 &mw->menu->zd->graph2d.copy);
  }
  dimN_R  =(dim1 *) malloc((graph.anzN_R+graph.anzN_Rp)*sizeof(dim1));
  for (i=0;i<graph.anzN_R;i++){
    dimN_R[i][0]=graph.dimN_R[i][0];
  }
  for (j=0;j<graph.anzN_Rp;i++,j++){
    dimN_R[i][0]=graph.dimN_Rp[j][0];
  }
  dimN_RxR=(dim1 *) malloc((graph.anzN_RxR+graph.anzN_RxRp)*sizeof(dim1));
  for (i=0;i<graph.anzN_RxR;i++){
    dimN_RxR[i][0]=graph.dimN_RxR[i][0];
  }
  for (j=0;j<graph.anzN_RxRp;i++,j++){
    dimN_RxR[i][0]=graph.dimN_RxRp[j][0];
  }
  dimN_RxRxRxR=(dim1 *) malloc((graph.anzN_RxRxRxR)*sizeof(dim1));
  for (i=0;i<graph.anzN_RxRxRxR;i++){
    dimN_RxRxRxR[i][0]=graph.dimN_RxRxRxR[i][0];
  }


  if (mw->menu->zd->graph2d.scaling){
    mw->menu->zd->graph2d.xmin=0;
    mw->menu->zd->graph2d.xmax=0;
    mw->menu->zd->graph2d.ymin=0;
    mw->menu->zd->graph2d.ymax=0;
    scalegraph2d(gp,mw->menu->zd->graph2d.draw,
		 graph.anzN_RxR+graph.anzN_RxRp,
		 dimN_RxR,&mw->menu->zd->graph2d.xmin,
		 &mw->menu->zd->graph2d.xmax,&mw->menu->zd->graph2d.ymin,
		 &mw->menu->zd->graph2d.ymax);
    /*** Scale graph can not assume that the scale of xmin and xmax
	 and ymin and ymax is already defined. solve this tomorrow. */
    scalegraph1d(gp1,mw->menu->zd->graph2d.draw1,
		 graph.anzN_R+graph.anzN_Rp,
		 dimN_R,
		 &mw->menu->zd->graph2d.xmin,&mw->menu->zd->graph2d.xmax,
		 &mw->menu->zd->graph2d.ymin,&mw->menu->zd->graph2d.ymax);
    scalegraph_tubes(gp2,mw->menu->zd->graph2d.draw2,
		 graph.anzN_RxRxRxR,
		 dimN_RxRxRxR,&mw->menu->zd->graph2d.xmin,
		 &mw->menu->zd->graph2d.xmax,&mw->menu->zd->graph2d.ymin,
		 &mw->menu->zd->graph2d.ymax);
    AddBorder(&mw->menu->zd->graph2d.xmin,&mw->menu->zd->graph2d.xmax,
	      &mw->menu->zd->graph2d.ymin,&mw->menu->zd->graph2d.ymax);
  }
  myclear();

  draw2dgraph(xsize,ysize,gp1,dimN_R,mw->menu->zd->graph2d.draw1,
	      graph.anzN_R+graph.anzN_Rp,mw->menu->zd->graph2d.comments,
	      mw->menu->zd->graph2d.xmin,mw->menu->zd->graph2d.xmax,
	      mw->menu->zd->graph2d.logx,
	      mw->menu->zd->graph2d.ymin,mw->menu->zd->graph2d.ymax,
	      mw->menu->zd->graph2d.logy,0,
	      mygrdat1);
  draw2dgraph(xsize,ysize,gp,dimN_RxR,mw->menu->zd->graph2d.draw,
	      graph.anzN_RxR+graph.anzN_RxRp,mw->menu->zd->graph2d.comments,
	      mw->menu->zd->graph2d.xmin,mw->menu->zd->graph2d.xmax,
	      mw->menu->zd->graph2d.logx,
	      mw->menu->zd->graph2d.ymin,mw->menu->zd->graph2d.ymax,
	      mw->menu->zd->graph2d.logy,1,
	      mygrdat);
  draw2dtubes(xsize,ysize,gp2,dimN_RxRxRxR,mw->menu->zd->graph2d.draw2,
	      graph.anzN_RxRxRxR,mw->menu->zd->graph2d.comments,
	      mw->menu->zd->graph2d.xmin,mw->menu->zd->graph2d.xmax,
	      mw->menu->zd->graph2d.logx,
	      mw->menu->zd->graph2d.ymin,mw->menu->zd->graph2d.ymax,
	      mw->menu->zd->graph2d.logy,1,
	      graph.dN_RxRxRxR);
  mw->menu->zd->time.newdata=0;
  free(mygrdat);
  free(mygrdat1);
  /* in this version gp's point to arrays and they do not allocate any
     memory to them. 
  for (i=0; i< graph.anzN_RxR+graph.anzN_RxR; i++)
    if (mw->menu->zd->graph2d.draw[i]) free(gp[i]);
    */
  if (gp!=NULL) free(gp);
  if (dimN_R!=NULL) free(dimN_R);
  if (gp1!=NULL) free(gp1);
  if (dimN_RxR!=NULL) free(dimN_RxR);
  
  myshow();
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}


int mycompare(const void *a,const void *b)
{
  int i=((int *) a)[1];
  int j=((int *) b)[1];
  return (i<j)? 1: (j<i)? -1:0;
}

void drawcontour(double *cont,int *length,int xmin,int ymin,int xmax, int ymax,
		 int xofs,int yofs,int nxsize, int nysize, int col1,
		 int col2)
{
  XPoint *points;
  int i,j,k,l,*sorted;
  double *contp;


  if (length[0]>0)
    {
      sorted=(int *) malloc(length[0]*2*sizeof(int));
      for (i=0;i<length[0];i++) {
	sorted[2*i]=i;
	sorted[2*i+1]=length[2*(i+1)];
      }
      qsort(sorted,length[0],2*sizeof(int),mycompare);
      l=sorted[1];
      points=(XPoint*) malloc(l*sizeof(XPoint));
      if (points==NULL) fprintf(stderr,"Memory error in drawcontour.\n");
      else
	for (j=0;j<length[0];j++)
	  {
	    for (i=0,contp=cont; i<sorted[2*j]; i++) contp+=2*length[2*(i+1)];
	    for (k=0;k<sorted[2*j+1];k++)
	      {
		points[k].x=xofs+(contp[k*2]-xmin+1)*nxsize/(xmax-xmin+1);
		points[k].y=yofs+nysize-(contp[k*2+1]-ymin+1)*nysize
		  /(ymax-ymin+1);
	      }
	    if (length[2*sorted[2*j]+3]== -4){
	      mypolygon(col1,points,sorted[2*j+1]);
	      myline_polygon(schwarz(),points,sorted[2*j+1]);}
	    else if (length[2*sorted[2*j]+3]== 4){
	      mypolygon(col2,points,sorted[2*j+1]);
	      myline_polygon(schwarz(),points,sorted[2*j+1]);
	    }
	    else
	      myline_polygon(gruen(),points,sorted[2*j+1]);
	  }
      free(points);
      free(sorted);
    }
}

void drawcontourline(double *cont,int *length,int xmin,int ymin,
		     int xmax, int ymax,
		     int xofs,int yofs,int nxsize, int nysize,int color)
{
  XPoint *points;
  int i,j,k,l,*sorted;
  double *contp;


  if (length[0]>0)
    {
      sorted=(int *) malloc(length[0]*2*sizeof(int));
      for (i=0;i<length[0];i++) {
	sorted[2*i]=i;
	sorted[2*i+1]=length[2*(i+1)];
      }
      qsort(sorted,length[0],2*sizeof(int),mycompare);
      l=sorted[1];
      points=(XPoint*) malloc(l*sizeof(XPoint));
      if (points==NULL) fprintf(stderr,"Memory error in drawcontour.\n");
      else
	for (j=0;j<length[0];j++)
	  {
	    for (i=0,contp=cont; i<sorted[2*j]; i++) contp+=2*length[2*(i+1)];
	    for (k=0;k<sorted[2*j+1];k++)
	      {
		points[k].x=xofs+(contp[k*2]-xmin+1)*nxsize/(xmax-xmin+1);
		points[k].y=yofs+nysize-(contp[k*2+1]-ymin+1)*nysize
		  /(ymax-ymin+1);
	      }
	    /*
	    if (length[2*sorted[2*j]+3]== -4){
	      myline_polygon(schwarz(),points,sorted[2*j+1]);}
	    else if (length[2*sorted[2*j]+3]== 4){
	      myline_polygon(schwarz(),points,sorted[2*j+1]);
	    }
	    else
	    */
	    myline_polygon(color,points,sorted[2*j+1]);
	  }
      free(points);
      free(sorted);
    }
}

void Cutdata(double *gpo,int ox,int oy,double *gpn,int nx,int ny,int center,
	     double *xofs, double *yofs)
{
  int x,y,cutx=0,cuty=0,cmxI,cmyI,ymax,ymin;
  double sum,cmx,cmy,*gppo;

  gppo=gpo;
  for (x=0;x<ox;x++) if (gpo[x*oy]<0) cuty=oy/2;
  for (y=0;y<oy;y++) if (gpo[y]<0) cutx=ox/2;
  if (center){
    sum=0;cmx=0;cmy=0; /* */ ymax=0; ymin=oy;
    for (x=0;x<ox;x++)
      for (y=0;y<oy;y++,gppo++)
	if (*gppo<0){
	  /* */ if (y<ymin) ymin=y;
	  /* */ if (y>ymax) ymax=y;
	  sum+=*gppo;
	  cmx+=*gppo*((x+cutx)%ox);
	  cmy+=*(gppo)*((y+cuty)%oy);
	}
    /* */ printf(" ymax=%i  ymin=%i\n",ymax,ymin);
    if (sum==0){
      fprintf(stderr,"Error in neuzeichnen.c:Cutdata: center finds no negative densities\n");
      cmxI=ox/2;cmyI=ox/2;
      cmx=cmy=0;
    }
    else {
      cmx /=sum;
      cmy /=sum;
      *xofs=cmx;
      *yofs=cmy;
      cmxI =floor(cmx);
      cmx -=cmxI;
      cmxI-=cutx;
      cmyI =floor(cmy);
      cmy -=cmyI;
      cmyI-=cuty;
    }
  }
  else {
    cmxI=floor(*xofs);cmyI=floor(*yofs);
    cmx=*xofs-cmxI;cmy=*yofs-cmyI;
  }
  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++)
      gpn[x*ny+y]=
	+(1-cmx)*(1-cmy) *gpo[((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy]
	+   cmx *(1-cmy) *gpo[((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy]
	+(1-cmx)*   cmy  *gpo[((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy]
	+   cmx *   cmy  *gpo[((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy];
  
}

void Cutdata_ve(double *gpo,int ox,int oy,double *gpn,int nx,int ny,
	     double xofs, double yofs)
{
  int x,y,cmxI,cmyI;
  double cmx,cmy,*gppo;

  gppo=gpo;

  cmxI=floor(xofs);cmyI=floor(yofs);
  cmx=xofs-cmxI;cmy=yofs-cmyI;

  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++){
      gpn[(x*ny+y)*2]=
	+(1-cmx)*(1-cmy) 
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*2]
	+   cmx *(1-cmy) 
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*2]
	+(1-cmx)*   cmy  
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*2]
	+   cmx *   cmy  
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*2];
      gpn[(x*ny+y)*2+1]=
	+(1-cmx)*(1-cmy) 
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*2+1]
	+   cmx *(1-cmy) 
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*2+1]
	+(1-cmx)*   cmy  
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*2+1]
	+   cmx *   cmy  
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*2+1];
    }
}

void Cutdata_te(double *gpo,int ox,int oy,double *gpn,int nx,int ny,
	     double xofs, double yofs)
{
  int x,y,cmxI,cmyI;
  double cmx,cmy,*gppo;

  gppo=gpo;

  cmxI=floor(xofs);cmyI=floor(yofs);
  cmx=xofs-cmxI;cmy=yofs-cmyI;

  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++){
      gpn[(x*ny+y)*3]=
	+(1-cmx)*(1-cmy) 
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3]
	+   cmx *(1-cmy) 
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3]
	+(1-cmx)*   cmy  
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3]
	+   cmx *   cmy  
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3];
      gpn[(x*ny+y)*3+1]=
	+(1-cmx)*(1-cmy) 
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3+1]
	+   cmx *(1-cmy) 
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3+1]
	+(1-cmx)*   cmy  
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3+1]
	+   cmx *   cmy  
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3+1];
      gpn[(x*ny+y)*3+2]=
	+(1-cmx)*(1-cmy) 
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3+2]
	+   cmx *(1-cmy) 
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+cmyI-ny/2+oy)%oy)*3+2]
	+(1-cmx)*   cmy  
	*gpo[(((x+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3+2]
	+   cmx *   cmy  
	*gpo[(((x+1+cmxI-nx/2+ox)%ox)*oy+(y+1+cmyI-ny/2+oy)%oy)*3+2];
    }
}

void neu_zeichnen2(Display *myd,mywindow *mw,int PS,char *psname, 
		   double PSwidth, double PSheight, int PSletter, 
		   int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,nxsize=1,nysize=1,no[2],no1[2],xofs=0,yofs=0,
    xlofs=0,ylofs=0,xlsize=0,ylsize=0,
    i,*length,linedist,
    xmin,xmax,ymin,ymax,contexist=0,boundary=0;
  double scale,scaletest,*cont,zmax,zmin;
  char comment[400];
  int graphcount;
  XPoint screenp[4];
  double *gp=NULL,*gp1=NULL,*sp=NULL,*gptmp=NULL,**vp;

   
  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }

  myclear();
  comment[0]='\0';

  /* some memory management that might be better kept in the menu
     section. Unfortunately an apropriate menu type has not yet been 
     defined. So for the time beeing I will keep it here. */
  /*
  if (mw->menu->zd->nz2.nocuts>=mw->menu->zd->nz2.maxcuts){
    mw->menu->zd->nz2.maxcuts=2*mw->menu->zd->nz2.nocuts;
    mw->menu->zd->nz2.cut=(double *) 
      realloc(mw->menu->zd->nz2.cut,mw->menu->zd->nz2.maxcuts*sizeof(double));
  }
  This should now all be contained in the double_arrp structres in menu.
  Aug 2004.
  */

  /* I can't be bothered to implement that propperly now... */
  if (mw->menu->zd->nz2.comments) ysize*=0.95;
  if (mw->menu->zd->nz2.DensityLegend) {
    xlsize=0.1*xsize;
    ylsize=ysize;
    xlofs=xsize-xlsize;
    ylofs=yofs;
    if (!mw->menu->zd->nz2.BarInset)
      xsize-=xlsize;
  }
  for (graphcount=0;graphcount<graph.anzNxN_R+graph.anzNxN_Rp;graphcount++)
    if (mw->menu->zd->nz2.draw[graphcount] && !mw->menu->zd->nz2.drawtwodensity){
      contexist=1;
      gp=NULL;
      if (graphcount<graph.anzNxN_R)
	getgraph3d(graph,NxN_R,&gp,
		   &no[0],graphcount,mw->menu->zd->nz2.rough,
		   &mw->menu->zd->nz2.copy[0]);
      else getgraph3d(graph,NxN_Rp,&gp,
		   &no[0],graph.anzNxN_R-graphcount,mw->menu->zd->nz2.rough,
		   &mw->menu->zd->nz2.copy[0]);
      if ((no[0]==0)||(no[1]==0)) contexist=0;
      
      if (contexist){
	if (mw->menu->zd->nz2.magnify){
	  gptmp=(double *) malloc(mw->menu->zd->nz2.NewXsize*
				  mw->menu->zd->nz2.NewYsize*sizeof(double));
	  Cutdata(gp,no[0],no[1],gptmp,
		  mw->menu->zd->nz2.NewXsize,mw->menu->zd->nz2.NewYsize,
		  mw->menu->zd->nz2.DoCenter,
		  &mw->menu->zd->nz2.xofs,&mw->menu->zd->nz2.yofs);
	  gp=gptmp;
	  no[0]=mw->menu->zd->nz2.NewXsize;
	  no[1]=mw->menu->zd->nz2.NewYsize;
	}
      
	xmin=0; xmax=no[0]; ymin=0; ymax=no[1];
	
	nxsize=xsize; nysize=ysize;
	if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
	  nysize=xsize*(ymax-ymin)/(xmax-xmin);
	else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
	if (nxsize<xsize) xofs=(xsize-nxsize)/2;
	if (nysize<ysize) yofs=(ysize-nysize)/2;
	
	
	if (mw->menu->zd->nz2.adjustcuts&&contexist){
	  scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  for (i=0;i<mw->menu->zd->nz2.nocuts;i++) 
	    mw->menu->zd->nz2.cut[i]=
	      zmin+(i+1)*(zmax-zmin)/(mw->menu->zd->nz2.nocuts+1);
	}
	if (mw->menu->zd->nz2.adjustdensity&&contexist){
	  scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  mw->menu->zd->nz2.density_min=zmin;
	  mw->menu->zd->nz2.density_max=zmax;
	}

	if ((mw->menu->zd->nz2.koord==1)&&contexist){
	  Koordinatensystem2D(schwarz(),mw->menu->zd->nz2.koordxmin,mw->menu->zd->nz2.koordxmax,mw->menu->zd->nz2.koordlogx,mw->menu->zd->nz2.koordymin,mw->menu->zd->nz2.koordymax,mw->menu->zd->nz2.koordlogy,&xofs,&yofs,
		      &nxsize,&nysize,"X","Y","","");
	}
	  /* At this point we should also allow a coordinate system to
	     be drawn around the density field. The standard values for the 
	     coordinates should simply be the dimensions of the field, but
	     one should be able to alter these values to correspond to
	     physical values. It would be preferable if this data was
	     available in the data stucture of NxN_R already so that it could
	     be defined in DefineGraph().
	  */


	if ((mw->menu->zd->nz2.density==1)&&contexist){
	  densityfield(no[0],no[1],mw->menu->zd->nz2.density_min,
		       mw->menu->zd->nz2.density_max,
		       ColorRampStart(),ColorRampEnd(),
		       gp,xofs,yofs,nxsize,nysize);


	  if (mw->menu->zd->nz2.DensityLegend){
	    sp=(double *)
	      malloc((ColorRampEnd()-ColorRampStart()+1)*sizeof(double));
	    for (i=0;i<(ColorRampEnd()-ColorRampStart()+1);i++)
	      sp[i]=mw->menu->zd->nz2.density_min
		+(mw->menu->zd->nz2.density_max
		  -mw->menu->zd->nz2.density_min)*i/
		(ColorRampEnd()-ColorRampStart());
	    if (nxsize<xsize) xlofs-=0.5*(xsize-nxsize);
	    if (nysize<ysize) ylsize=nysize;
	    ylofs=yofs;
	    if (mw->menu->zd->nz2.BarInset)
	      xlofs=xofs+mw->menu->zd->nz2.BarXOfs*xsize;

	    densityfield(1,ColorRampEnd()-ColorRampStart()+1,
			 mw->menu->zd->nz2.density_min,
			 mw->menu->zd->nz2.density_max,
			 ColorRampStart(),ColorRampEnd(),
			 sp,xlofs-
			 0.5*mw->menu->zd->nz2.BarWidth*xlsize,
			 ylofs+mw->menu->zd->nz2.BarYOfs*ylsize,
			 2*mw->menu->zd->nz2.BarWidth*xlsize,
			 mw->menu->zd->nz2.BarHeight*ylsize);
	    free(sp);
	    CoordLineY(schwarz(), mw->menu->zd->nz2.density_min,
		       mw->menu->zd->nz2.density_max,0 /*log ?*/,
		       xlofs+mw->menu->zd->nz2.BarWidth*xlsize,
		       ylofs+mw->menu->zd->nz2.BarYOfs*ylsize,
		       (1-mw->menu->zd->nz2.BarWidth)*xlsize,
		       mw->menu->zd->nz2.BarHeight*ylsize);
	  }
	  if (UserComment==NULL)
	    sprintf(comment,"Min=%g, Max=%g, ",mw->menu->zd->nz2.density_min,
		    mw->menu->zd->nz2.density_max);
	}
	for (i=0;(i<mw->menu->zd->nz2.nocuts)&&(contexist!=0);i++)
	  {
	    contour(gp,no[0],no[1],
		    mw->menu->zd->nz2.cut[i],&cont,&length,0,NULL,
		    &boundary);
	    switch (mw->menu->zd->nz2.density)
	      {
	      case 0:
		screenp[0].x=xofs+nxsize/(no[0]+1);
		screenp[0].y=yofs+nysize/(no[1]+1);
		screenp[1].x=xofs+nxsize/(no[0]+1);
		screenp[1].y=yofs+nysize-nysize/(no[1]+1);
		screenp[2].x=xofs+nxsize-nxsize/(no[0]+1);
		screenp[2].y=yofs+nysize-nysize/(no[1]+1);
		screenp[3].x=xofs+nxsize-nxsize/(no[0]+1);
		screenp[3].y=yofs+nysize/(no[1]+1);
		if (boundary==0 && gp[0]>mw->menu->zd->nz2.cut[i])
		  mypolygon(blau(),screenp,4);
		else
		  mypolygon(gelb(),screenp,4);

		/* Not quite right! depends on the color of the edge states!!! */
		drawcontour(cont,length,xmin,ymin,
			    xmax,ymax,
			    xofs,yofs,nxsize,nysize,gelb(),blau());
		break;
	      case 1 :
	      case 2 :
		drawcontourline(cont,length,xmin,ymin,
				xmax,ymax,
				xofs,yofs,nxsize,nysize,schwarz());
		break;
	      default: 
		fprintf(stderr,"error in neuzeichnen.c:neu_zeichnen2() density=%i",
			mw->menu->zd->nz2.density);
		
	      }
	    free(cont);
	    free(length);
	  }
	if (mw->menu->zd->nz2.magnify) free(gptmp);
      }
      /*  free(gp); gp=NULL; Not needed again! */
  }

  /* Sept 2010: new data display of two simultanous densities, colored according to a color field. This needs to be finished. Define new part in menu structure. Need a menu to select the two graphs to compare, then draw the graphs here!*/
  if (mw->menu->zd->nz2.drawtwodensity){
    contexist=1;
    gp=NULL;
    for (i=0;i<graph.anzNxN_R+graph.anzNxN_Rp;i++){
      if (mw->menu->zd->nz2.draw[i]) {
	mw->menu->zd->nz2.d1=i;
	i++;
	break;
      }
    }
    for (;i<graph.anzNxN_R+graph.anzNxN_Rp;i++){
      if (mw->menu->zd->nz2.draw[i]){
	mw->menu->zd->nz2.d2=i;
	break;
      }
    }
    if (i==graph.anzNxN_R+graph.anzNxN_Rp) mw->menu->zd->nz2.d2=mw->menu->zd->nz2.d2=0;
    if (mw->menu->zd->nz2.d1<graph.anzNxN_R)
      getgraph3d(graph,NxN_R,&gp,
		 &no[0],mw->menu->zd->nz2.d1,mw->menu->zd->nz2.rough,
		 &mw->menu->zd->nz2.copy[0]);
    else getgraph3d(graph,NxN_Rp,&gp,
		    &no[0],mw->menu->zd->nz2.d1-graph.anzNxN_R,mw->menu->zd->nz2.rough,
		    &mw->menu->zd->nz2.copy[0]);
    if (mw->menu->zd->nz2.d2<graph.anzNxN_R)
      getgraph3d(graph,NxN_R,&gp1,
		 &no1[0],mw->menu->zd->nz2.d2,mw->menu->zd->nz2.rough,
		 &mw->menu->zd->nz2.copy[0]);
    else getgraph3d(graph,NxN_Rp,&gp1,
		    &no1[0],mw->menu->zd->nz2.d2-graph.anzNxN_R,mw->menu->zd->nz2.rough,
		    &mw->menu->zd->nz2.copy[0]);
    if ((no[0]!=no1[0])||(no[1]!=no1[1])) {
      contexist=0;
      printf("Error in drawing twodensity graph: graphs have different dimensions!\n x: %i=%i, y: %i=%i\n",no[0],no1[0],no[1],no1[1]);
    }
    
    if (contexist){
	if (mw->menu->zd->nz2.magnify){
	  gptmp=(double *) malloc(mw->menu->zd->nz2.NewXsize*
				  mw->menu->zd->nz2.NewYsize*sizeof(double));
	  Cutdata(gp,no[0],no[1],gptmp,
		  mw->menu->zd->nz2.NewXsize,mw->menu->zd->nz2.NewYsize,
		  mw->menu->zd->nz2.DoCenter,
		  &mw->menu->zd->nz2.xofs,&mw->menu->zd->nz2.yofs);
	  gp=gptmp;
	  gptmp=(double *) malloc(mw->menu->zd->nz2.NewXsize*
				  mw->menu->zd->nz2.NewYsize*sizeof(double));
	  Cutdata(gp1,no[0],no[1],gptmp,
		  mw->menu->zd->nz2.NewXsize,mw->menu->zd->nz2.NewYsize,
		  mw->menu->zd->nz2.DoCenter,
		  &mw->menu->zd->nz2.xofs,&mw->menu->zd->nz2.yofs);
	  gp1=gptmp;
	  no[0]=mw->menu->zd->nz2.NewXsize;
	  no[1]=mw->menu->zd->nz2.NewYsize;
	}
      
	xmin=0; xmax=no[0]; ymin=0; ymax=no[1];
	
	nxsize=xsize; nysize=ysize;
	if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
	  nysize=xsize*(ymax-ymin)/(xmax-xmin);
	else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
	if (nxsize<xsize) xofs=(xsize-nxsize)/2;
	if (nysize<ysize) yofs=(ysize-nysize)/2;
	
	
	if (mw->menu->zd->nz2.adjustcuts&&contexist){
	  scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  for (i=0;i<mw->menu->zd->nz2.nocuts;i++) 
	    mw->menu->zd->nz2.cut[i]=
	      zmin+(i+1)*(zmax-zmin)/(mw->menu->zd->nz2.nocuts+1);
	  scaleScalar(gp1,no[0]*no[1],&zmin,&zmax);
	  for (i=0;i<mw->menu->zd->nz2.nocuts1;i++) 
	    mw->menu->zd->nz2.cut1[i]=
	      zmin+(i+1)*(zmax-zmin)/(mw->menu->zd->nz2.nocuts+1);
	}
	if (mw->menu->zd->nz2.adjustdensity&&contexist){
	  scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  mw->menu->zd->nz2.density1_min=zmin;
	  mw->menu->zd->nz2.density1_max=zmax;
	  scaleScalar(gp1,no[0]*no[1],&zmin,&zmax);
	  mw->menu->zd->nz2.density2_min=zmin;
	  mw->menu->zd->nz2.density2_max=zmax;
	}

	if ((mw->menu->zd->nz2.koord==1)&&contexist){
	  Koordinatensystem2D(schwarz(),mw->menu->zd->nz2.koordxmin,mw->menu->zd->nz2.koordxmax,mw->menu->zd->nz2.koordlogx,mw->menu->zd->nz2.koordymin,mw->menu->zd->nz2.koordymax,mw->menu->zd->nz2.koordlogy,&xofs,&yofs,
		      &nxsize,&nysize,"X","Y","","");
	}
	  /* At this point we should also allow a coordinate system to
	     be drawn around the density field. The standard values for the 
	     coordinates should simply be the dimensions of the field, but
	     one should be able to alter these values to correspond to
	     physical values. It would be preferable if this data was
	     available in the data stucture of NxN_R already so that it could
	     be defined in DefineGraph().
	  */


	if ((mw->menu->zd->nz2.drawtwodensity==1)&&contexist){
	  twodensityfield(no[0],no[1],
			  mw->menu->zd->nz2.density1_min,
			  mw->menu->zd->nz2.density1_max,
			  mw->menu->zd->nz2.density2_min,
			  mw->menu->zd->nz2.density2_max,
			  ColorFieldStart(),ColorFieldX(),ColorFieldY(),
			  gp,gp1,xofs,yofs,nxsize,nysize);

	  if (UserComment==NULL)
	    sprintf(comment,"Min1=%g, Max1=%g, Min2=%g, Max2=%g, ",
		    mw->menu->zd->nz2.density1_min,
		    mw->menu->zd->nz2.density1_max,
		    mw->menu->zd->nz2.density2_min,
		    mw->menu->zd->nz2.density2_max);
	}
	for (i=0;(i<mw->menu->zd->nz2.nocuts)&&(contexist!=0);i++)
	  {
	    contour(gp,no[0],no[1],
		    mw->menu->zd->nz2.cut[i],&cont,&length,0,NULL,
		    &boundary);
	    drawcontourline(cont,length,xmin,ymin,
			    xmax,ymax,
			    xofs,yofs,nxsize,nysize,gruen());
	    free(cont);
	    free(length);
	  }
	for (i=0;(i<mw->menu->zd->nz2.nocuts1)&&(contexist!=0);i++)
	  {
	    contour(gp1,no[0],no[1],
		    mw->menu->zd->nz2.cut1[i],&cont,&length,0,NULL,
		    &boundary);
	    drawcontourline(cont,length,xmin,ymin,
			    xmax,ymax,
			    xofs,yofs,nxsize,nysize,rot());
	    free(cont);
	    free(length);
	  }
	if (mw->menu->zd->nz2.magnify) free(gptmp);
    } 
  }
  
  scale=0;
  contexist=0;
  vp=(double **) malloc(graph.anzNxN_RxR*sizeof(double *));
  for (i=0;i<graph.anzNxN_RxR;i++){
    vp[i]=NULL;
    if (mw->menu->zd->nz2.draw_ufield[i])
      {
	contexist=1;
	getgraph3d(graph,NxN_RxR,&(vp[i]),
		   &no[0],i,mw->menu->zd->nz2.rough,
		   &mw->menu->zd->nz2.copy[0]);
	if (mw->menu->zd->nz2.magnify){
	  gptmp=(double *) malloc(2*mw->menu->zd->nz2.NewXsize*
				  mw->menu->zd->nz2.NewYsize*sizeof(double));
	  Cutdata_ve(vp[i],no[0],no[1],gptmp,
		  mw->menu->zd->nz2.NewXsize,mw->menu->zd->nz2.NewYsize,
		  mw->menu->zd->nz2.xofs,mw->menu->zd->nz2.yofs);
	  vp[i]=gptmp;
	  no[0]=mw->menu->zd->nz2.NewXsize;
	  no[1]=mw->menu->zd->nz2.NewYsize;

	}

	xmin=0; xmax=no[0]; ymin=0; ymax=no[1];

	nxsize=xsize; nysize=ysize;
	if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
	  nysize=xsize*(ymax-ymin)/(xmax-xmin);
	else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
	if (nxsize<xsize) xofs=(xsize-nxsize)/2;
	if (nysize<ysize) yofs=(ysize-nysize)/2;

	scaletest=vec_scale(no[0]*no[1],vp[i]);
	if (scaletest>scale) scale=scaletest;
      }
  }
  if ((mw->menu->zd->nz2.koord==1)&&contexist){
    Koordinatensystem2D(schwarz(),mw->menu->zd->nz2.koordxmin,mw->menu->zd->nz2.koordxmax,mw->menu->zd->nz2.koordlogx,mw->menu->zd->nz2.koordymin,mw->menu->zd->nz2.koordymax,mw->menu->zd->nz2.koordlogy,&xofs,&yofs,
			&nxsize,&nysize,"X","Y","","");
  }

  if (contexist!=0){
    if (UserComment==0) 
      sprintf(comment,"%s u max = %g",comment,scale);
    scale/=nxsize/(xmax-xmin+1);
    for (i=0;i<graph.anzNxN_RxR;i++){
      if (mw->menu->zd->nz2.draw_ufield[i])
	{
	  ufield(no[1],xmin,xmax,
		 ymin,ymax,
		 vp[i],scale,
		 xofs,yofs,nxsize,nysize);
	  if (mw->menu->zd->nz2.magnify) free(vp[i]);/* allocated as gptmp */
	}
    }
  }
  free(vp);

  vp=(double **) malloc(graph.anzNxN_RxRxRt*sizeof(double *));
  scale=0;
  contexist=0;
  for (i=0;i<graph.anzNxN_RxRxRt;i++){
    vp[i]=NULL;
    if (mw->menu->zd->nz2.draw_tfield[i])
      {
	contexist=1;
	getgraph3d(graph,NxN_RxRxRt,&(vp[i]),
		   &no[0],i,mw->menu->zd->nz2.rough,
		   &mw->menu->zd->nz2.copy[0]);
	if (mw->menu->zd->nz2.magnify){
	  gptmp=(double *) malloc(3*mw->menu->zd->nz2.NewXsize*
				  mw->menu->zd->nz2.NewYsize*sizeof(double));
	  Cutdata_te(vp[i],no[0],no[1],gptmp,
		  mw->menu->zd->nz2.NewXsize,mw->menu->zd->nz2.NewYsize,
		  mw->menu->zd->nz2.xofs,mw->menu->zd->nz2.yofs);
	  vp[i]=gptmp;
	  no[0]=mw->menu->zd->nz2.NewXsize;
	  no[1]=mw->menu->zd->nz2.NewYsize;

	}

	xmin=0; xmax=no[0]; ymin=0; ymax=no[1];

	nxsize=xsize; nysize=ysize;
	if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
	  nysize=xsize*(ymax-ymin)/(xmax-xmin);
	else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
	if (nxsize<xsize) xofs=(xsize-nxsize)/2;
	if (nysize<ysize) yofs=(ysize-nysize)/2;

	scaletest=te_scale(no[0]*no[1],vp[i]);
	if (scaletest>scale) scale=scaletest;
      }
  }
  if (contexist!=0){
    if (UserComment==NULL)
      sprintf(comment,"%s t max = %g",comment,scale);
    scale/=nxsize/(xmax-xmin+1);
    for (i=0;i<graph.anzNxN_RxRxRt;i++){
      if (mw->menu->zd->nz2.draw_tfield[i])
	{
	  tfield(no[1],xmin,xmax,
		 ymin,ymax,
		 vp[i],scale,
		 xofs,yofs,nxsize,nysize,mw->menu->zd->nz2.rough_t[i]);
	  if (mw->menu->zd->nz2.magnify) free(vp[i]);/* allocated as gptmp */
	}
    }
  }
  free(vp);

  if (mw->menu->zd->nz2.comments)
    {
      if (UserComment!=NULL) sprintf(comment,"%s",UserComment);
      myselectfont(0,"My nice comment",1.0*xsize/4);
      linedist=1.0*ysize/30;
      mytext(schwarz(),xsize/2,yofs+nysize+linedist,comment,1);
    }
  myshow();
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}

void zeichnen_3d(Display *myd,mywindow *mw,int PS,char *psname, 
		   double PSwidth, double PSheight, int PSletter, 
		   int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,nxsize=1,nysize=1,no3d[3],no[2],xofs=0,yofs=0,
    xlofs=0,ylofs=0,xlsize=0,ylsize=0,
    i,*length,linedist,
    xmin,xmax,ymin,ymax,contexist=0,boundary=0;
  double scale,scaletest,*cont,zmax,zmin;
  char comment[400];
  int graphcount;
  XPoint screenp[4];
  double *gp=NULL,*sp=NULL,*gptmp=NULL,**vp;

   
  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }

  myclear();
  comment[0]='\0';

  /* some memory management that might be better kept in the menu
     section. Unfortunately an apropriate menu type has not yet been 
     defined. So for the time beeing I will keep it here. */
  /*
  if (mw->menu->zd->nz2.nocuts>=mw->menu->zd->nz2.maxcuts){
    mw->menu->zd->nz2.maxcuts=2*mw->menu->zd->nz2.nocuts;
    mw->menu->zd->nz2.cut=(double *) 
      realloc(mw->menu->zd->nz2.cut,mw->menu->zd->nz2.maxcuts*sizeof(double));
  }
  This should now all be contained in the double_arrp structres in menu.
  Aug 2004.
  */

  /* I can't be bothered to implement that properly now... */
  if (mw->menu->zd->p3d.comments) ysize*=0.95;
  for (graphcount=0;graphcount<graph.anzNxNxN_R;graphcount++)
    if (mw->menu->zd->p3d.drawX[graphcount]){
      contexist=1;
      gp=NULL;
      getgraph3d(graph,NxNxN_R,&gp,
		 &no3d[0],graphcount,1,
		     &mw->menu->zd->p3d.copy[0]);
      no[0]=no3d[1];
      no[1]=no3d[2];
      if ((no[0]==0)||(no[1]==0)) contexist=0;
      
      if (contexist){
	if (mw->menu->zd->p3d.magnify){
	  gptmp=(double *) malloc(mw->menu->zd->p3d.NewXsize*
				  mw->menu->zd->p3d.NewYsize*sizeof(double));
	  Cutdata(gp,no[0],no[1],gptmp,
		  mw->menu->zd->p3d.NewXsize,mw->menu->zd->p3d.NewYsize,
		  0 /* center */,
		  &mw->menu->zd->p3d.xofs,&mw->menu->zd->p3d.yofs);
	  gp=gptmp;
	  no[0]=mw->menu->zd->p3d.NewXsize;
	  no[1]=mw->menu->zd->p3d.NewYsize;
	}
      
	xmin=0; xmax=no[0]; ymin=0; ymax=no[1];
	
	nxsize=xsize; nysize=ysize;
	if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
	  nysize=xsize*(ymax-ymin)/(xmax-xmin);
	else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
	if (nxsize<xsize) xofs=(xsize-nxsize)/2;
	if (nysize<ysize) yofs=(ysize-nysize)/2;
	
	
	if (mw->menu->zd->p3d.adjustdensity&&contexist){
	  scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  mw->menu->zd->p3d.density_min=zmin;
	  mw->menu->zd->p3d.density_max=zmax;
	}
	if ((mw->menu->zd->p3d.density==1)&&contexist){
	  densityfield(no[0],no[1],mw->menu->zd->p3d.density_min,
		       mw->menu->zd->p3d.density_max,
		       ColorRampStart(),ColorRampEnd(),
		       gp,xofs,yofs,nxsize,nysize);

	  /* At this point we should also allow a coordinate system to
	     be drawn around the density field. The standard values for the 
	     coordinates should simply be the dimensions of the field, but
	     one should be able to alter these values to correspond to
	     physical values. It would be preferable if this data was
	     available in the data stucture of NxN_R already so that it could
	     be defined in DefineGraph().
	  */
 
	  if (UserComment==NULL)
	    sprintf(comment,"Min=%g, Max=%g, ",mw->menu->zd->p3d.density_min,
		    mw->menu->zd->p3d.density_max);
	}
	if (mw->menu->zd->p3d.magnify) free(gptmp);
      }
      /*  free(gp); gp=NULL; Not needed again! */
  }
  
  
  if (mw->menu->zd->p3d.comments)
    {
      if (UserComment!=NULL) sprintf(comment,"%s",UserComment);
      myselectfont(0,"My nice comment",1.0*xsize/4);
      linedist=1.0*ysize/30;
      mytext(schwarz(),xsize/2,yofs+nysize+linedist,comment,1);
    }
  myshow();
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}

int Contour_print(char *psname, double PSwidth, double PSheight, 
                   int PSletter, int PSlandscape,
		   int nocuts,  double *cut, double *gp, int nox, int noy,
		   int comments, const char *comment, int magnify, 
		   int center,
		   int NewXsize, int NewYsize,double xofs, double yofs,
		   int adjustcuts,
		   int adjustdensity, double density_min, 
		   double density_max, int density)

{
  int xsize,ysize,nxsize=1,nysize=1,no[2],xofsw=0,yofsw=0,
    i,*length,linedist,boundary=0,
    xmin,xmax,ymin,ymax,contexist=1;
  double *cont,zmax,zmin;
  XPoint screenp[4];
  double *gptmp=NULL;

  no[0]=nox;
  no[1]=noy;

  if (initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize))
    return 1;
 
  /*setPS(); myclear(); should be unnecessary. */

  if ((no[0]==0)||(no[1]==0)) contexist=0;

  if (comments) ysize*=0.95;
      
  if (contexist){
    if (magnify){
      gptmp=(double *) malloc(NewXsize*NewYsize*sizeof(double));
      Cutdata(gp,no[0],no[1],gptmp,NewXsize,NewYsize,center,&xofs,&yofs);
      gp=gptmp;
      no[0]=NewXsize;
      no[1]=NewYsize;
    }
    
    xmin=0; xmax=no[0]; ymin=0; ymax=no[1];
	
    nxsize=xsize; nysize=ysize;
    if (ysize>xsize*(ymax-ymin)/(xmax-xmin))
      nysize=xsize*(ymax-ymin)/(xmax-xmin);
    else nxsize=ysize*(xmax-xmin)/(ymax-ymin);
    if (nxsize<xsize) xofsw=(xsize-nxsize)/2;
    if (nysize<ysize) yofsw=(ysize-nysize)/2;
    
	
    if (adjustcuts&&contexist){
      scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
      for (i=0;i<nocuts;i++) 
	cut[i]=
	  zmin+(i+1)*(zmax-zmin)/(nocuts+1);
    }
    if (adjustdensity&&contexist){
      scaleScalar(gp,no[0]*no[1],&zmin,&zmax);
	  density_min=zmin;
	  density_max=zmax;
    }
    if ((density==1)&&contexist){
      densityfield(no[0],no[1],density_min,
		   density_max,
		   ColorRampStart(),ColorRampEnd(),
		   gp,xofs,yofs,nxsize,nysize);
    }
    for (i=0;(i<nocuts)&&(contexist!=0);i++)
      {
	contour(gp,no[0],no[1],
		cut[i],&cont,&length,0,NULL,&boundary);
	switch (density)
	  {
	  case 0:
	    screenp[0].x=xofsw+nxsize/(no[0]+1);
	    screenp[0].y=yofsw+nysize/(no[1]+1);
	    screenp[1].x=xofsw+nxsize/(no[0]+1);
	    screenp[1].y=yofsw+nysize-nysize/(no[1]+1);
	    screenp[2].x=xofsw+nxsize-nxsize/(no[0]+1);
	    screenp[2].y=yofsw+nysize-nysize/(no[1]+1);
	    screenp[3].x=xofsw+nxsize-nxsize/(no[0]+1);
	    screenp[3].y=yofsw+nysize/(no[1]+1);
	    if (boundary==0 && gp[0]>cut[i])
	      mypolygon(blau(),screenp,4);
	    else
	      mypolygon(gelb(),screenp,4);

	    /* Not quite right! depends on the color of the edge states!!! */
	    drawcontour(cont,length,xmin,ymin,
			xmax,ymax,
			xofs,yofs,nxsize,nysize,gelb(),blau());
	    break;
	  case 1 :
	  case 2 :
	    drawcontourline(cont,length,xmin,ymin,
			    xmax,ymax,
			    xofs,yofs,nxsize,nysize,schwarz());
	    break;
	  default: 
	    fprintf(stderr,"error in neuzeichnen.c:neu_zeichnen2() density=%i",
		    density);
	    
	  }
	free(cont);
	free(length);
      }
    if (magnify) free(gptmp);
  }
  

  if (comments)
    {
      myselectfont(0,"My nice comment",1.0*xsize/4);
      linedist=1.0*ysize/30;
      mytext(schwarz(),xsize/2,yofs+nysize+linedist,comment,1);
    }
  myshow();
  setX();
  return 0;
}

void neu_zeichnen_contour3d(Display *myd,mywindow *mw,int PS,char *psname, 
			    double PSwidth, double PSheight, int PSletter, 
			    int PSlandscape, grstr graph)
{
  int i,xsize,ysize,no[3],xmin,xmax,ymin,ymax,zmax,zmin;
  Contour **c=NULL;
  int nocont=0;
  
  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }
  if (graph.anzcont3d>0)
    c=(Contour **) malloc(graph.anzcont3d*sizeof(Contour *));
  for (i=0;i<graph.anzcont3d;i++)
    if (mw->menu->zd->nzc.getgraph[i]){
      getgraph3d(graph,cont3d,(double **)&(c[nocont]),
		 &no[0],i,mw->menu->zd->nzc.rough,
		 &mw->menu->zd->nzc.copy[0]);
      nocont++;
    }

  xmin=0; xmax=no[0]*mw->menu->zd->nz1.rough; 
  ymin=0; ymax=no[1]*mw->menu->zd->nz1.rough;
  zmin=0; zmax=no[2]*mw->menu->zd->nz1.rough;

  
  draw3dcontour(xsize,ysize,c,nocont,
	      mw->menu->zd->nzc.newdata,
	      mw->menu->zd->nzc.persp_factor,mw->menu->zd->nzc.scalefakt, 
	      mw->menu->zd->nzc.xangle,mw->menu->zd->nzc.yangle,0,
	      xmax,xmin,ymax,ymin,zmin,zmax);

  if (c!=NULL) free(c);
  mw->menu->zd->nz1.newdata=0;
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}



#ifdef notdef
/* this should be rewritten using fftw!*/
void neu_zeichnen_grfft(Display *myd,mywindow *mw,int PS,char *psname, double PSwidth, double PSheight, int PSletter, int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,no[2],xmin,xmax,ymin,ymax,i,x,y,m,ct,tvar;
  double *field;
  double average;
  static int initialized=0;
  static double **gp;
  
  if (initialized==0){
    gp=(double **) malloc(graph.anzNxN_R*sizeof(double *));
    for (i=0;i<graph.anzNxN_R;i++) gp[i]=NULL;
    initialized=1;
  }
  if (mw->menu->zd->nz1.newdata)
    {
      for (i=0;i<graph.anzNxN_R;i++)
	if (mw->menu->zd->nz1.getgraph[i]!=0)
	  /* adapt to new standard !
	  getgraph3d(graph,NxN_R,&(gp[i]),
		     &no[0],i,mw->menu->zd->nz1.rough,
		     mw->menu->zd->nz1.copyx,mw->menu->zd->nz1.copyy);
	else
	  if(gp[i]!=NULL) {free(gp[i]); gp[i]=NULL;}
	  */

      /* Checking that size=2^n */
      tvar=no[0];
      for (m=1,ct=0;m<16;ct+=tvar&1,tvar>>=1,m++);
      if (ct!=1) return;
      tvar=no[1];
      for (m=1,ct=0;m<16;ct+=tvar&1,tvar>>=1,m++);
      if (ct!=1) return;
      
#ifndef notdef
      field=malloc(no[0]*no[1]*2*sizeof(double));
      for (i=0;i<graph.anzNxN_R;i++)
	{
	  if (mw->menu->zd->nz1.getgraph[i])
	    {
	      average=0;
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++)
		  average+=gp[i][x*no[1]+y];
	      average/=no[0]*no[1];
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++){
		  field[2*(x*no[1]+y)]=gp[i][x*no[1]+y]-average;
		  field[2*(x*no[1]+y)+1]=0;
		}
	      fourn(field-1,no-1,2,1);
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++)
		  gp[i][x*no[1]+y]=
		    pow2(field[2*(((x+no[0]/2)%no[0])*no[1]+
				      (y+no[1]/2)%no[1])])
			 +pow2(field[2*(((x+no[0]/2)%no[0])*no[1]+
				       (y+no[1]/2)%no[1])+1]);

	    }
	}
      free(field);

#else
      field=malloc(no[0]*no[1]*4*sizeof(double));
      for (i=0;i<graph.anzS;i++)
	{
	  if (mw->menu->zd->nz1.getgraph[i])
	    {
	      average=0;
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++)
		  average+=
		    gp[i][x*no[1]+y];
	      average/=no[0]*no[1];
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++){
		  field[2*(x*no[1]*2+y)]=field[2*(x*no[1]*2+2*no[1]-1-y)]=
		    gp[i][x*no[1]+y]-average;
		  field[2*(x*no[1]*2+y)+1]=field[2*(x*no[1]*2+2*no[1]-1-y)+1]=0;
		}
	      no[1]*=2;
	      fourn(field-1,no-1,2,1);
	      for (x=0;x<no[0]*no[1];x++) {
		field[2*x]=pow2(field[2*x])+pow2(field[2*x+1]);
		field[2*x+1]=0;
	      }
 	      fourn(field-1,no-1,2,-1);
	      no[1]/=2;
	      for (x=0;x<no[0];x++)
		for (y=0;y<no[1]; y++)
		  gp[i][x*no[1]+y]=
				field[2*(((x+no[0]/2)%no[0])*no[1]*2+
					 ((y*2+no[1])%(no[1]*2)))]+
		                field[2*(((x+no[0]/2)%no[0])*no[1]*2+
				         ((y*2+1+no[1])%(no[1]*2)))];
	    }
	}
      free(field);
#endif	  
      if (mw->menu->zd->nz1.Autoscale||mw->menu->zd->nz1.scaleonce)
	scalegraph(gp,mw->menu->zd->nz1.getgraph,graph.anzNxN_R,
		   no[0]*no[1],
		   &mw->menu->zd->nz1.zmin,&mw->menu->zd->nz1.zmax);
      mw->menu->zd->nz1.scaleonce=0;
    }
  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }
  xmin=0; xmax=no[0]*mw->menu->zd->nz1.rough; ymin=0; ymax=no[1]*mw->menu->zd->nz1.rough;
  draw3dgraph(xsize,ysize,graph,gp,mw->menu->zd->nz1.getgraph,
	      no[0],no[1],
	      &mw->menu->zd->nz1.NEWDIMENSIONS,mw->menu->zd->nz1.newdata,mw->menu->zd->nz1.comments,
	      /*get_time()*/0,
	      mw->menu->zd->nz1.persp_factor,mw->menu->zd->nz1.scalefakt, 
	      mw->menu->zd->nz1.xangle,mw->menu->zd->nz1.yangle,0,
	      xmax,xmin,ymax,ymin,
	      mw->menu->zd->nz1.zmin,mw->menu->zd->nz1.zmax,mw->menu->zd->nz1.zfact,
	      &(mw->menu->zd->nz1.me),&(mw->menu->zd->nz1.meinit));
  mw->menu->zd->nz1.newdata=0;
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}

void neu_zeichnen_Structfact(Display *myd,mywindow *mw,int PS,char *psname, double PSwidth, double PSheight, int PSletter, int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,i,Slen,nograph;
  double **gp=NULL,ymin,ymax,ymintmp,ymaxtmp;

  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }
  nograph=0;
  for (i=0;i<graph.anzNxN_R;i++)
    {
      if (mw->menu->zd->nz1.getgraph[i])
	{
	  nograph++;
	  gp=(double **) realloc(gp,nograph*sizeof(double *));
	  gp[nograph-1]=get_structfact(graph,i,&Slen);
	  get_minmax(gp[nograph-1],Slen,&ymintmp,&ymaxtmp);
	  if (nograph==1) {ymin=ymintmp;ymax=ymaxtmp;}
	  else{
	    if (ymin>ymintmp) ymin=ymintmp;
	    if (ymax<ymaxtmp) ymax=ymaxtmp;
	  }
	}
    }
  if (nograph==0) return;
  /* adapt to new standerd.
  draw2dgraph(xsize,ysize,gp,Slen,nograph,mw->menu->zd->nz1.comments,
	      Slen-1,0,ymax,ymin,0);
	      */
  mw->menu->zd->nz1.newdata=0;
  for (i=0;i<nograph;i++) free(gp[i]);
  free(gp);

  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}

void neu_zeichnen_data(Display *myd,mywindow *mw,int PS,char *psname, double PSwidth, double PSheight, int PSletter, int PSlandscape,
		   grstr graph)
{
  int xsize,ysize,linedist,line=1,Slen,size;
  char comment[200],l1s[30],l0s[30];
  double d0,d2[3],l0,l1,tet,tmp,S_0,S_1,S_2,S_3,myl,kav,d1[2];


  if (PS)
    {
      initPS(psname,PSwidth,PSheight,PSletter,PSlandscape,&xsize,&ysize);
      setPS();
    }
  else
    {
      setXwindow(myd,mw->zeichenw,mw->mypixmap,mw->zeichengc,
	    mw->doublebuffer,mw->zeichenmaxx,mw->zeichenmaxy);
      xsize=mw->zeichenmaxx;
      ysize=mw->zeichenmaxy;
      mw->newdraw=0;
    }
  myclear();
  myselectfont(0,"My nice comment",1.0*xsize/3);
  linedist=1.0*ysize/13;
/*  Here we need to check that we can get*/
  get_fftdata(graph,0,&d0,&kav,d1,d2); 
  l0=(d2[0]+d2[2])/2
    -sqrt(pow2((d2[0]+d2[2])/2)-d2[0]*d2[2]
	  +pow2(d2[1]));
  l1=(d2[0]+d2[2])/2
    +sqrt(pow2((d2[0]+d2[2])/2)-d2[0]*d2[2]
	  +pow2(d2[1]));
  if ((l0-l1)*(d2[0]-d2[2])<0) 
    { tmp=l0;l0=l1;l1=tmp;}
  if (d2[1]==0) tet=0;
  else tet=asin(1/sqrt(pow2((d2[2]-l0))+pow2(d2[1]))*d2[1]);

  sprintf(comment,"d0=%e",d0);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"d2=(%e,%e,%e)",d2[0],d2[1],d2[2]);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  if (l0==0) sprintf(l0s,"infinite"); else sprintf(l0s,"%e",sqrt(d0/l0));
  if (l1==0) sprintf(l1s,"infinite"); else sprintf(l1s,"%e",sqrt(d0/l1));
  sprintf(comment,"l1=%s l2=%s tet=%e",l0s,l1s,tet/M_PI*180);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  /**/
  get_derivdata(graph,3,d2); 
  d0=get_av2(graph,0);

  l0=(d2[0]+d2[2])/2
    -sqrt(pow2((d2[0]+d2[2])/2)-d2[0]*d2[2]
	  +pow2(d2[1]));
  l1=(d2[0]+d2[2])/2
    +sqrt(pow2((d2[0]+d2[2])/2)-d2[0]*d2[2]
	  +pow2(d2[1]));
  if ((l0-l1)*(d2[0]-d2[2])<0) 
    { tmp=l0;l0=l1;l1=tmp;}
  if (d2[1]==0) tet=0;
  else tet=asin(1/sqrt(pow2((d2[2]-l0))+pow2(d2[1]))*d2[1]);

  mytext(schwarz(),xsize/2,line*linedist,"Calculated from derivatives:",1);line++;
  sprintf(comment,"d0=%e",d0);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"d2=(%e,%e,%e)",d2[0],d2[1],d2[2]);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  myl=sqrt(d0/(l0+l1));
  if (l0==0) sprintf(l0s,"infinite"); 
  else sprintf(l0s,"%e",sqrt(d0/2./l0));
  if (l1==0) sprintf(l1s,"infinite"); 
  else sprintf(l1s,"%e",sqrt(d0/2./l1));
  sprintf(comment,"l=%e l1=%s l2=%s tet=%e",myl,
	  l0s,l1s,tet/M_PI*180);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  
  get_Struct_data(graph,0,&S_0,&S_1,&S_2,&S_3,&Slen,&size);
  sprintf(comment,"S_0=%e",S_0);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"S_1=%e L*S_0/S_1=%e",S_1,size*S_0/S_1);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"S_2=%e L(S_0/S_2)^(1/2)=%e L*S_1/S_2=%e ",S_2,
	  size*sqrt(S_0/S_2),size*S_1/S_2);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"S_3=%e L(S_0/S_3)^(1/3)=%e L(S_1/S_3)^(1/2)=%e",S_3,
	  size*pow(S_0/S_3,1./3),size*sqrt(S_1/S_3) );
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  sprintf(comment,"L(S_1/S_3)^(1/2)/l=%e",
	  size*sqrt(S_1/S_3)/myl);
  mytext(schwarz(),xsize/2,line*linedist,comment,1);line++;
  

  myshow();
  if (PS)
    {
      setX();
    }
  else XFlush(myd);
}

#endif







