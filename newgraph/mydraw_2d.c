/***************************************************************************\
graph library
Copyright (C) 1992-2003  Alexander Wagner and Johannes Schlosser

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


#include "mydraw_2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

static Display   
    *mydisplay=NULL;

static Drawable  
    mydrawable,
    mydrawabledouble, 
    mydrawto;

static GC    
    mygc;


static int       
    shape = Nonconvex,
    mode  = CoordModeOrigin,
    myXactive = 1,
    myXmode,
    xsize = 1000,
    ysize = 1000,
    PSxsize,
    PSysize,
    Xxsize,
    Xysize;

#define  PSscale 0.056
#define  maxPSxsizeA4 10000.0
#define  maxPSysizeA4 (maxPSxsizeA4*29.5/21) /* DIN A4 */

#define  maxPSxsizeLetter 10278.0
#define  maxPSysizeLetter (maxPSxsizeLetter*11/8.5) /* American Letter */

static int maxPSxsize,maxPSysize,PSlandscape;
static int PSxmin,PSymin,PSxmax,PSymax; 
static FILE  *PSfile;
static char  *PSfilename;
static fpos_t BBoxPos;

/****   Color-Map spezifisches (Benutzung privater Color-Cells) ****/


#define PIXELS 1024
                        
static  unsigned long           mypixels [PIXELS];

	struct PScolor{ double red, green, blue;};
/*	typedef struct PScolor[PIXELS] PSmap;*/
	
static 	struct PScolor myPSmap[PIXELS];

int _schwarz=0,_weiss=0,_rot=0,_gruen=0,_blau=0,_gelb=0;


void (*mydrawline)(int color, int x1, int y1, int x2,int y2);
void (*mypolygon)(int color, XPoint *points, int n);
void (*mypolygon_line)(int color, int color_line, XPoint *points, int n);
void (*myline_polygon)(int color_line, XPoint *points, int n);
void (*mytriangle_2)(int color, int color_line, XPoint *points);
void (*mycircle)(int color, int x, int y, int r);
void (*myfilledcircle)(int color, int x, int y, int r);
void (*densityfield)(int nox, int noy,double zmin, double zmax,
		     int  colmin, int colmax, double *f,
		     int xofs, int yofs, int xsize, int ysize);

void (*mytext)(int color,int x, int y,const char *text, int orient);
void (*myselectfont)(int font, char *testtext, double length);
void (*mysetfontdirection)(double alpha);

void (*myclear)();
void (*myshow)();




int schwarz() {return _schwarz;}
int weiss(){return _weiss;}
int rot(){return _rot;}
int gruen(){return _gruen;}
int blau(){return _blau;}
int gelb(){return _gelb;}
void set_schwarz(int f) {_schwarz=f;}
void set_weiss(int f){_weiss=f;}
void set_rot(int f){_rot=f;}
void set_gruen(int f){_gruen=f;}
void set_blau(int f){_blau=f;}
void set_gelb(int f){_gelb=f;}


void PSxlimit(int x)
{
  if (x<PSxmin) PSxmin=x;
  if (x>PSxmax) PSxmax=x;
}
void PSylimit(int y)
{
  if (y<PSymin) PSymin=y;
  if (y>PSymax) PSymax=y;
}

void Visual_Information(int *Class,unsigned long *red_mask,
			unsigned long *green_mask, unsigned long *blue_mask,
			Visual **visual)
{
  XVisualInfo visualinfotemp, *visualinfo;
  Visual *def_visual;
  int no;

  visualinfo =
  XGetVisualInfo(mydisplay,VisualNoMask /*| VisualIDMask | VisualScreenMask |
		 VisualDepthMask | VisualClassMask | VisualRedMaskMask |
		 VisualGreenMaskMask | VisualBlueMaskMask | 
		 VisualColormapSizeMask | VisualBitsPerRGBMask*/ 
		 /*VisualAllMask*/, &visualinfotemp, &no);
  if (visualinfo==NULL) {
    printf(" NO VISUALS FOUND \n");
    return;
  }
  /*  printf("%i visuals found\n",no);
  printf("GrayScale=%i StaticGray=%i PseudoColor=%i StaticColor=%i "
	 "DirectColor=%i TrueColor=%i\n",GrayScale,StaticGray,PseudoColor,
	   StaticColor,DirectColor,TrueColor);
  for (i=0; i<no; i++){
    printf("Screen: %i\n",visualinfo[i].screen); 
    printf("Depth: %u\n",visualinfo[i].depth); 
    printf("Class: %i\n",visualinfo[i].class);
    printf("red_mask: %ul\n",visualinfo[i].red_mask); 
    printf("green_mask: %ul\n",visualinfo[i].green_mask); 
    printf("blue_mask: %ul\n",visualinfo[i].blue_mask); 
    printf("Colormap Size: %i\n",visualinfo[i].colormap_size); 
    printf("Bits per RGB: %i\n",visualinfo[i].bits_per_rgb);
    } */
  def_visual=DefaultVisual(mydisplay,0);
  /*  printf("The DefaultVisual is:\n");
  printf("Class: %i\n",def_visual->class);
  printf("red_mask: %ul\n",def_visual->red_mask); 
  printf("green_mask: %ul\n",def_visual->green_mask); 
  printf("blue_mask: %ul\n",def_visual->blue_mask); 
  printf("Bits per RGB: %i\n",def_visual->bits_per_rgb);*/
  *Class=def_visual->class;
  *red_mask=def_visual->red_mask;
  *green_mask=def_visual->green_mask;
  *blue_mask=def_visual->blue_mask;
  *visual=def_visual;
}

static  Bool   ColorInitialized=0;
static  Colormap                mycmap;
static int NumberOfCol,RampStart=0,RampEnd=0;
static  int Class;

int ColorRampStart(){ return RampStart;}
int ColorRampEnd() {return RampEnd;}
int NumberOfColors(){return NumberOfCol;}
void SetColorRampStart(int i){RampStart=i;}
void SetColorRampEnd  (int i){RampEnd =i;}


int GetFreeColors(){
#define PLANES 8 
  static  Bool                    mycontig=0;
  static  unsigned long           myplane_masks [PLANES];
  static 	unsigned int            mynplanes       = 0;
  static  Status                  myresult;

  static  unsigned long red_mask,green_mask,blue_mask;
  static  Visual *myvisual;

  Visual_Information(&Class,&red_mask,&green_mask,&blue_mask,&myvisual);
  mycontig	= 1;
  mycmap	= DefaultColormap (mydisplay,0);
  if (Class==4){
    NumberOfCol=PIXELS;
    ColorInitialized = 1;
    return NumberOfCol;
  }
  else {
    for (NumberOfCol=256;NumberOfCol>0;NumberOfCol--){
      myresult= XAllocColorCells(mydisplay, mycmap, mycontig, myplane_masks,
				 mynplanes, mypixels, NumberOfCol);
      if (myresult!=0) {
	ColorInitialized = 1;
	return NumberOfCol;
      }
    }
  }
  printf("Allocation of colors failed in GetFreeColors()\n");
  NumberOfCol=0;
  return NumberOfCol;
}

void setcolormap(XColor *map, int n)
{
  int i;

  if (mydisplay!=NULL){
    if (!ColorInitialized){
      GetFreeColors();
    }
    if (n<=NumberOfCol) {
      for (i=0; i<n; i++) {
	myPSmap[i].red	= map [i].red  /65535.0;
	myPSmap[i].green= map [i].green/65535.0;
	myPSmap[i].blue = map [i].blue /65535.0;
	if (Class !=4){
	  map [i].pixel	= mypixels [i];
	  map [i].flags	= DoRed | DoGreen | DoBlue;
	}
	else {
	  if (!XAllocColor(mydisplay,mycmap,&(map[i])))
	    printf("True Color allocation error");
	  mypixels[i]=map[i].pixel;
	}
      }
      if (Class !=4)
	XStoreColors(mydisplay, mycmap, map, n);
    }
  }
  else { /* for PostScript usage only */
    for (i=0; i<n; i++) {
      myPSmap[i].red  = 1.525902e-05*map [i].red;
      myPSmap[i].green= 1.525902e-05*map [i].green;
      myPSmap[i].blue = 1.525902e-05*map [i].blue;
    }
  }
}

void myclearPS(){
  rewind(PSfile);
  fprintf(PSfile,"%%!PS-Adobe-3.0 EPSF-3.0\n"); /* Hier Header erzeugen */
  fprintf(PSfile,"%%%%Creator: automatically produced by mydraw.c"
	  " (author Johannes Schlosser/Alexander Wagner)\n");
  fgetpos(PSfile,&BBoxPos);
  if (PSlandscape)
    fprintf(PSfile,"%%%%BoundingBox: %07.3f %07.3f %07.3f %07.3f             ",
	    0.0,PSscale*(maxPSxsize-xsize),PSscale*ysize,PSscale*maxPSxsize);
  else
    fprintf(PSfile,"%%%%BoundingBox: %07.3f %07.3f %07.3f %07.3f             ",
	    0.0,0.0,PSscale*xsize,PSscale*maxPSysize);
  fprintf(PSfile,"\n/Pol {lineto lineto setrgbcolor fill} def\n");
  fprintf(PSfile,"\n/Tgl {newpath moveto lineto lineto setrgbcolor "
            "gsave fill grestore 0.0 0.0 0.0 setrgbcolor stroke} def\n");
  fprintf(PSfile,"\n/P4 {newpath moveto lineto lineto lineto setrgbcolor "
            "gsave fill grestore setrgbcolor stroke} def\n");
  fprintf(PSfile,"gsave\n");
  if (PSlandscape)
    fprintf(PSfile,"90 rotate\n");
  else
    fprintf(PSfile,"0 %07.3f translate\n",PSscale*maxPSysize);
  fprintf(PSfile,"%7.3f %7.3f scale\n", PSscale, PSscale);
  fprintf(PSfile,"\n");
  fprintf(PSfile,	"/rightshow\n"
	  "{ dup stringwidth pop\n"
	  "  neg 0 rmoveto\n"
	  "  show\n"
	  "} def\n\n"
	  "/centershow\n"
	  "{ dup stringwidth pop\n"
	  "  -2 div 0 rmoveto\n"
	  "  show\n"
	  "} def\n\n");
  
  fprintf(PSfile,"5 setlinewidth\n");
  PSxmin=PSxsize;PSymin=PSysize;PSxmax=0;PSymax=0;
}

void myclearX(){
  XSetForeground(mydisplay, mygc, mypixels[0]);
  XFillRectangle(mydisplay, mydrawto, mygc, 0, 0, xsize, ysize);
}

void myshowX(){
  if (myXmode == doublebufferX)
    XCopyArea(mydisplay, mydrawto, mydrawable, mygc, 0, 0,
	      xsize, ysize, 0, 0);    
}

void myshowPS(){
  int i;
  fprintf(PSfile,"showpage\n");
  fsetpos(PSfile,&BBoxPos);
  if (PSlandscape)
    fprintf(PSfile,"%%%%BoundingBox: %07.3f %07.3f %07.3f %07.3f             ",
	    0.0,PSscale*(maxPSxsize-xsize),PSscale*ysize,PSscale*maxPSxsize);
  else
    fprintf(PSfile,"%%%%BoundingBox: %07.3f %07.3f %07.3f %07.3f             ",
	    PSscale*PSxmin,PSscale*(maxPSysize-PSymax),
	    PSscale*PSxmax,PSscale*(maxPSysize-PSymin));
  i=fclose(PSfile);
  if (i!=0) printf("Could not close PostScript file! error %i\n",i);
  PSfile=NULL;
}                                             

void mydrawlineX(int color, int x1, int y1, int x2,int y2){
  XSetForeground(mydisplay, mygc, mypixels[color]);
  XDrawLine(mydisplay, mydrawto, mygc, x1, y1, x2, y2);
}

void mypolygonX(int color, XPoint *points, int n){
	XSetForeground(mydisplay, mygc, mypixels[color]);
	XFillPolygon(mydisplay, mydrawto, mygc, points, n, shape, mode);
}

void mypolygon_lineX(int color, int color_line, XPoint *points, int n){
	XSetForeground(mydisplay, mygc, mypixels[color]);
	XFillPolygon(mydisplay, mydrawto, mygc, points, n, shape, mode);
	XSetForeground(mydisplay,mygc,mypixels[color_line]);
	XDrawLines(mydisplay, mydrawto, mygc, points, n, mode);
}

void myline_polygonX(int color_line, XPoint *points, int n){
	XSetForeground(mydisplay,mygc,mypixels[color_line]);
	XDrawLines(mydisplay, mydrawto, mygc, points, n, mode);
}

void mycircleX(int color, int x, int y, int r){
        XSetForeground(mydisplay, mygc, mypixels[color]);
	XDrawArc(mydisplay,mydrawto, mygc, x-r, y-r, 2*r, 2*r, 0, 23040);
      }

void myfilledcircleX(int color, int x, int y, int r){
        XSetForeground(mydisplay, mygc, mypixels[color]);
	XFillArc(mydisplay,mydrawto, mygc, x-r, y-r, 2*r, 2*r, 0, 23040);
      }

void mydrawlinePS(int color, int x1, int y1, int x2,int y2){
  PSxlimit(x1);
  PSxlimit(x2);
  PSylimit(y1);
  PSylimit(y2);
  fprintf(PSfile,"newpath\n");
  fprintf(PSfile,"%d %d moveto\n",x1,-y1);
  fprintf(PSfile,"%d %d lineto\n",x2,-y2);
  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue);	
  fprintf(PSfile,"stroke\n");
}

void mycirclePS(int color, int x, int y, int r){
  PSxlimit(x+r);
  PSxlimit(x-r);
  PSylimit(y+r);
  PSylimit(y-r);
  fprintf(PSfile,"newpath\n");
  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue);
  fprintf(PSfile,"%d %d %d 0 360 arc\n", x, -y, r);
  fprintf(PSfile,"stroke\n");
}

void myfilledcirclePS(int color, int x, int y, int r){
  PSxlimit(x+r);
  PSxlimit(x-r);
  PSylimit(y+r);
  PSylimit(y-r);
	fprintf(PSfile,"newpath\n");
	fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
		myPSmap[color].red,
		myPSmap[color].green,
		myPSmap[color].blue);
	fprintf(PSfile,"%d %d %d 0 360 arc\n", x, -y, r);
        fprintf(PSfile,"fill\n");
}



void mypolygonPS(int color, XPoint *points, int n){
  int i;

  for (i=0;i<n;i++){PSxlimit(points[i].x);PSylimit(points[i].y);}
  fprintf(PSfile,"newpath\n");
  fprintf(PSfile,"%d %d moveto\n",points[0].x,-points[0].y);
  for (i=1;i<n;i++)	
    fprintf(PSfile,"%d %d lineto\n",points[i].x,-points[i].y);	
  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue);	
  fprintf(PSfile,"fill\n");	
}

void mypolygon_line4PS(int color, int color_line, XPoint *points){
  int i;
  for (i=0;i<4;i++){PSxlimit(points[i].x);PSylimit(points[i].y);}
  fprintf(PSfile,"%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f"
          " %d %d %d %d %d %d %d %d  P4\n",
	  myPSmap[color_line].red,
	  myPSmap[color_line].green,
	  myPSmap[color_line].blue,
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue,
	  points[0].x, -points[0].y,
	  points[1].x, -points[1].y,
	  points[2].x, -points[2].y,
	  points[3].x, -points[3].y);
}


void mypolygon_linePS(int color, int color_line, XPoint *points, int n){
	int i;
	for (i=0;i<n;i++){PSxlimit(points[i].x);PSylimit(points[i].y);}
	if (n==4)
	  mypolygon_line4PS(color,color_line,points);
	else{
	  fprintf(PSfile,"newpath\n");
	  fprintf(PSfile,"%d %d moveto\n",points[0].x,-points[0].y);
	  for (i=1;i<n;i++)	
	    fprintf(PSfile,"%d %d lineto\n",points[i].x,-points[i].y);
	  fprintf(PSfile,"closepath\n");
	  fprintf(PSfile,"gsave\n");  	
	  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
		  myPSmap[color].red,
		  myPSmap[color].green,
		  myPSmap[color].blue);	
	  fprintf(PSfile,"fill\n");
	  fprintf(PSfile,"grestore\n");
	  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
		  myPSmap[color_line].red,
		  myPSmap[color_line].green,
		  myPSmap[color_line].blue);	
	  fprintf(PSfile,"stroke\n");
	}
}

void myline_polygonPS(int color_line, XPoint *points, int n){
  int i;
  for (i=0;i<n;i++){PSxlimit(points[i].x);PSylimit(points[i].y);}
  fprintf(PSfile,"newpath\n");
  fprintf(PSfile,"%d %d moveto\n",points[0].x,-points[0].y);
  for (i=1;i<n;i++)	
    fprintf(PSfile,"%d %d lineto\n",points[i].x,-points[i].y);
/*  fprintf(PSfile,"closepath\n");*/
  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
	  myPSmap[color_line].red,
	  myPSmap[color_line].green,
	  myPSmap[color_line].blue);	
  fprintf(PSfile,"stroke\n");
}

void mytriangle_2PS(int color, int color_line, XPoint *points){
  int i;
  for (i=0;i<3;i++){PSxlimit(points[i].x);PSylimit(points[i].y);}
  fprintf(PSfile,"%5.3f %5.3f %5.3f %d %d %d %d %d %d  Tgl\n",
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue,
	  points[1].x, -points[1].y,
	  points[0].x, -points[0].y,
	  points[2].x, -points[2].y);
}

void mytriangle_2X(int color, int color_line, XPoint *points){
  mypolygon(color, points, 3);
  mydrawline(color_line, points[0].x, points[0].y, points[1].x, points[1].y);
  mydrawline(color_line, points[0].x, points[0].y, points[2].x, points[2].y);
}

void densityfieldX(int nox, int noy,
	    double zmin, double zmax,int  colmin, int colmax, double *f,
	    int xofs, int yofs, int xsize, int ysize)
{
  int i,j,x,x1,y,y1;
  XPoint xp[4];
  double d,dx,dy;
  
  if (zmin!=zmax)
    {
      yofs+=ysize;
      dx=1.*xsize/(nox+1);
      dy=1.*ysize/(noy+1);
      for (i=0;i<nox;i++){
	x=xofs+dx*(i+0.5);
	x1=xofs+dx*(i+1.5);
	for (j=0;j<noy;j++)
	  {
	    y= yofs-dy*(j+0.5);
	    y1=yofs-dy*(j+1.5);
	    xp[0].x=x;
	    xp[0].y=y;
	    xp[1].x=x;
	    xp[1].y=y1;
	    xp[2].x=x1;
	    xp[2].y=y1;
	    xp[3].x=x1;
	    xp[3].y=y;
	    d=(f[i*noy+j]-zmin)/(zmax-zmin);
	    if (isnormal(d)){
	    if (d<0) d=0;
	    else if (d>1) d=1;
	    } else d=0;
	    mypolygon(colmin+d*(colmax-colmin),xp,4);
	  }
      }
    }
}

void densityfieldPS(int nox, int noy,
	    double zmin, double zmax,int  colmin, int colmax, double *f,
	    int xofs, int yofs, int xsize, int ysize)
{
  int i,j,col;
  int r,g,b;
  double d;

  if (zmin!=zmax)
    {
      yofs+=ysize;
      PSxlimit(xofs+0.5*xsize/(nox+1));
      PSylimit(yofs-0.5*ysize/(noy+1));
      PSxlimit(xofs+xsize*(1-.5/(nox+1)));
      PSylimit(yofs-ysize*(1-.5/(noy+1)));
	       
      fprintf(PSfile,"gsave\n");
      fprintf(PSfile,"/picstr %i string def\n",nox*3);
      fprintf(PSfile,"%e %e translate\n",
	      xofs+0.5*xsize/(nox+1),-yofs+0.5*ysize/(noy+1));
      fprintf(PSfile,"%e %e scale\n",
	      xsize*(1-1.0/(nox+1)),ysize*(1-1.0/(noy+1)));
      fprintf(PSfile,"%% Image geometry\n %i %i 8\n",nox,noy);
      fprintf(PSfile,"%% Transformation Matrix\n [%i 0 0 %i 0 0]\n",
	      nox,noy);
      fprintf(PSfile,"{currentfile picstr readhexstring pop}\n");
      fprintf(PSfile,"false 3\ncolorimage\n");
      for (j=0;j<noy;j++){
	for (i=0;i<nox;i++){
	    d=(f[i*noy+j]-zmin)/(zmax-zmin);
	    if (d<0) d=0;
	    else if (d>1) d=1;
	    col=colmin+d*(colmax-colmin);
	    r=myPSmap[col].red*255;
	    g=myPSmap[col].green*255;
	    b=myPSmap[col].blue*255;
	    fprintf(PSfile,"%02x%02x%02x",r,g,b);
	}
	fprintf(PSfile,"\n");
      }
      fprintf(PSfile,"grestore\n\n");
    }
}


struct			font_struct {	int l;
				XFontStruct *fontstr;
			    };
#define NOFONTS 6
static int FOUNDFONTS=0;
struct font_struct	fonts[NOFONTS];
XFontStruct             *aktfont;
char			*fontname[NOFONTS] = {
					"5x7",
					"5x8",
					"6x9",
					"6x10",
					"8x13",
					"9x15",
					};

char                    *PSfontname[]={
                                       "Times-Roman"
				     };
double                  psfontalpha;

void myXfontsinit(void)
{
  
 char **allfonts;
  char *myfontname[5] = {"5x8","6x9","6x10","8x13","9x15"};
  char fontbasename[] = {"*"};
  char fontpattern[300];

  static int initialized=0;
  int i,j;
  const char *teststr="X";

  if (initialized==0){
    initialized=1;

    sprintf(fontpattern, "-*-%s-medium-r-*", fontbasename);
    allfonts = XListFonts(mydisplay, fontpattern, NOFONTS, &FOUNDFONTS);

    for (i=0,j=0; i < FOUNDFONTS; i++)
      { /* Bekannte Fonts : fixed ,9x15, 8x13,6x10*/ 
	fonts[j].fontstr = XLoadQueryFont ( mydisplay,allfonts[i]);
	if (fonts[j].fontstr == 0) 
	  printf("Der Font %s konnte nicht gefunden werden. \n",fontname[i]);
	else {
	  fonts[j].l = XTextWidth(fonts[i].fontstr,teststr,strlen(teststr));
	  j++;
	}
      }
    FOUNDFONTS=j;
  }
}

void mytext_X(int color,int x, int y,const char *text, int orient)
{
  int length;
  
  XSetForeground(mydisplay, mygc, mypixels[color]);
  length = XTextWidth( aktfont, text, strlen(text));
  if (orient == 1)
    XDrawString( mydisplay, mydrawto, mygc, x-length/2,y, text, strlen(text));
  else if(orient == 2)
    XDrawString( mydisplay, mydrawto, mygc, x-length,y, text, strlen(text));
  else XDrawString( mydisplay, mydrawto, mygc, x,y, text, strlen(text));
}

void myselectfont_X(int font, char *testtext, double length)
{
  int i,lengthX;

  lengthX = length/strlen(testtext);
  for (i=0; (lengthX > fonts[i].l) && (i<FOUNDFONTS); i++);
  if ( i == 0 ) i=1;
  XSetFont (mydisplay, mygc, fonts[i-1].fontstr->fid);
  aktfont = fonts[i-1].fontstr;
}

void mysetfontdirection_X(double alpha)
{/* Es gibt keine gedrehten Koordinaten unter X... */
}


void mytext_PS(int color,int x, int y,const char *text, int orient)
{
  PSxlimit(x);PSylimit(y);
/* NOT REALLY ACCEPTABLE.... */
  fprintf(PSfile,"gsave\n");
  fprintf(PSfile,"%5.3f %5.3f %5.3f setrgbcolor\n",
	  myPSmap[color].red,
	  myPSmap[color].green,
	  myPSmap[color].blue);	

  fprintf(PSfile," %d %d moveto\n",x,-y);
  fprintf(PSfile," %f rotate\n",psfontalpha);
    switch(orient)
      {
      case 2: /* to act mark */
	fprintf(PSfile,"(");
	fprintf(PSfile,text);
	fprintf(PSfile,") rightshow\n");
	break;
	
      case 1: /* centered */
	fprintf(PSfile,"(");
	fprintf(PSfile,text);
	fprintf(PSfile,") centershow\n");
	break;
	
      case 0: /* from act. mark */
	fprintf(PSfile,"(");
	fprintf(PSfile,text);
	fprintf(PSfile,") show\n");
	break;
      }
  fprintf(PSfile,"grestore\n");
}

void myselectfont_PS(int font, char *testtext, double length)
{
  fprintf(PSfile,"/");
  fprintf(PSfile,PSfontname[font]);
  fprintf(PSfile,	" findfont dup\n"
	  "30 scalefont\n"
	  "setfont\n"
	  "(");
  fprintf(PSfile,	testtext);
  fprintf(PSfile,	") stringwidth pop\n"
	  "30 %10.5f  mul exch div scalefont\n"
	  "setfont\n", length);
  
}

void mysetfontdirection_PS(double alpha)
{
  psfontalpha=alpha;
}

void setdisplay(Display *md){
  mydisplay=md;
}

void setXwindow(Display *md, Drawable mydr, Drawable mydoubledr, GC myGC, 
		int mode, int _xsize, int _ysize){
  Xxsize=_xsize;
  Xysize=_ysize;
  xsize=Xxsize;
  ysize=Xysize;
  setdisplay(md);
  mydrawable=mydr;
  mydrawabledouble=mydoubledr;
  mydrawto=mydrawable;
  if (mode == doublebufferX) mydrawto=mydrawabledouble;
  mygc=myGC;
  myXmode=mode;
}

void initX(Display *md, Drawable mydr, Drawable mydoubledr, GC myGC, 
	   int mode, int _xsize, int _ysize){
	setXwindow(md,mydr,mydoubledr,myGC,mode,_xsize,_ysize);
	myXfontsinit();
	setX();
}


void setPS(){
    xsize=PSxsize;
    ysize=PSysize;
    mydrawline = mydrawlinePS;
    mypolygon = mypolygonPS;
    mypolygon_line = mypolygon_linePS;
    myline_polygon = myline_polygonPS;
    mytriangle_2 = mytriangle_2PS;
    mycircle = mycirclePS;
    myfilledcircle = myfilledcirclePS;
    densityfield = densityfieldPS;
    mytext = mytext_PS;
    myselectfont = myselectfont_PS;
    mysetfontdirection = mysetfontdirection_PS;
    myclear = myclearPS;
    myshow = myshowPS;
    myXactive = 0;
    myclear();
}

int initPS(char *filename, double _xsize, double _ysize, 
	    int PSletter, int PSlandscape_,int *_PSxsize, int *_PSysize){
  switch (PSletter){
  case 0:
    if (PSlandscape_){
      PSysize=maxPSysize=_xsize*maxPSxsizeA4;
      PSxsize=maxPSxsize=_ysize*maxPSysizeA4;
    }
    else{
      PSxsize=maxPSxsize=_xsize*maxPSxsizeA4;
      PSysize=maxPSysize=_ysize*maxPSysizeA4;
    }
    break;
  case 1:
    if (PSlandscape_){
      PSysize=maxPSysize=_xsize*maxPSxsizeLetter;
      PSxsize=maxPSxsize=_ysize*maxPSysizeLetter;
    }
    else{
      PSxsize=maxPSxsize=_xsize*maxPSxsizeLetter;
      PSysize=maxPSysize=_ysize*maxPSysizeLetter;
    }
    break;
  default: 
    fprintf(stderr,"mydrawc:initPS:Error format %i not 0=A4 or 1=letter\n",
	    PSletter);
    PSxsize=_xsize*maxPSxsizeA4;
    PSysize=_ysize*maxPSysizeA4;
    break;
  }
  PSlandscape=PSlandscape_;
  PSfilename=filename;
  PSfile=fopen(PSfilename,"w");
  if (PSfile==NULL) {
    fprintf(stderr,"Could not open file %s! Waiting... ",PSfilename);
    sleep(4);
    PSfile=fopen(PSfilename,"w");
    if (PSfile==NULL) {
      fprintf(stderr," failed again!\n");
      return 1;
    }
    else fprintf(stderr," O.K. now.\n");
  }
    
  setPS();
  *_PSxsize=PSxsize;
  *_PSysize=PSysize;
  return 0;
}

void setX(){
    xsize=Xxsize;
    ysize=Xysize;
    mydrawline = mydrawlineX;
    mypolygon = mypolygonX;
    mypolygon_line = mypolygon_lineX;
    myline_polygon = myline_polygonX;
    mytriangle_2 = mytriangle_2X;
    mycircle = mycircleX;
    myfilledcircle = myfilledcircleX;
    densityfield = densityfieldX;
    mytext = mytext_X;
    myselectfont = myselectfont_X;
    mysetfontdirection = mysetfontdirection_X;
    myclear = myclearX;
    myshow = myshowX;
    myXactive = 1;
}



