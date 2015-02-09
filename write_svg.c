/*----------------------------------------------------------------------------

  ELSD - Ellipse and Line Segment Detector

  Copyright (c) 2012 viorica patraucean (vpatrauc@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "elsd.h"
#include "write_svg.h"
#include "process_curve.h"
/*----------------------------------------------------------------------------*/
FILE* init_svg(char * filename, unsigned int xsize, unsigned int ysize )
{
  FILE * svg;

  /* open file */
  svg = fopen(filename,"w");
  if( svg == NULL ) error("Error: unable to open SVG output file.");

  /* write SVG header */
  fprintf(svg,"<?xml version=\"1.0\" standalone=\"no\"?>\n");
  fprintf(svg,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
  fprintf(svg," \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf(svg,"<svg width=\"%upx\" height=\"%upx\" ",xsize,ysize);
  fprintf(svg,"version=\"1.1\"\n xmlns=\"http://www.w3.org/2000/svg\" ");
  fprintf(svg,"xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");
  //fprintf(svg,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"",1,1,xsize,ysize);
  //fprintf(svg, " fill=\"none\" stroke=\"black\" stroke-width=\"%d\" />\n",1);

return svg;
 
}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void fclose_svg(FILE *svg)
{
 /* close SVG file */
  fprintf(svg,"</svg>\n");
  if( fclose(svg) == EOF )
    error("Error: unable to close file while writing SVG file.");
}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void write_svg_line(FILE *svg, double *lin,int smooth)
{
  /* write line segment */

  if (smooth)
    {
     lin[0] = lin[0]*1.25;
     lin[1] = lin[1]*1.25;
     lin[2] = lin[2]*1.25;
     lin[3] = lin[3]*1.25;
    }
  fprintf(svg,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" ",
              lin[0],lin[1],lin[2],lin[3]);
  fprintf(svg," fill=\"none\" stroke =\"green\" stroke-width=\"%d\" />\n",1);
  //fprintf(svg," fill=\"none\" stroke =\"black\" stroke-width=\"%d\" />\n",1);
}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void write_svg_circle(FILE *fe, FILE *svg,double *param,int *pext,int smooth)
{
  int fa=0,fs=1;
  double ang_start,ang_end;
  
  if (smooth) 
    {
      param[0] *= 1.25; param[1] *= 1.25;
      param[2] *= 1.25; param[3] *= 1.25;
      pext[0] *= 1.25; pext[1] *= 1.25;
      pext[2] *= 1.25; pext[3] *= 1.25;
    }
  
  ang_start = atan2(pext[1]-param[1],pext[0]-param[0]); 
  ang_end = atan2(pext[3]-param[1],pext[2]-param[0]);

  double C = M_2__PI * param[2];
  if (angle_diff(ang_start,ang_end)<M_2__PI*SQRT2/C && angle_diff_signed(ang_start,ang_end)>0) 
    {
      fprintf(svg,"<ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\" stroke=\"blue\" stroke-width=\"%d\" fill=\"none\" /> \n",
              param[0],param[1], param[2], param[3],1); 
    }
  else 
    {
      double x1,y1,x2,y2;
      x1 = param[2]*cos(ang_start) + param[0];
      y1 = param[2]*sin(ang_start) + param[1];
      x2 = param[2]*cos(ang_end) + param[0];
      y2 = param[2]*sin(ang_end) + param[1];
      if (ang_start<0) ang_start+=M_2__PI;
      if (ang_end<0) ang_end+=M_2__PI;
  
      if (ang_end<ang_start) ang_end += M_2__PI;

      if ((ang_end-ang_start)>M_PI) fa = 1;
      fprintf(svg,"<path d=\"M %f,%f A%f,%f %f %d,%d %f,%f\"",
	x1,y1,param[2],param[2],0.0,fa,fs,x2,y2);
      fprintf(svg," fill=\"none\" stroke =\"blue\" stroke-width=\"%d\" />\n",1); 
      //fprintf(svg," fill=\"none\" stroke =\"black\" stroke-width=\"%d\" />\n",1);   
    }
  fprintf(fe,"%f %f %f %f %f %f %f\n",param[0],param[1],param[2],param[3],param[4],ang_start,ang_end);

}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void write_svg_ellipse(FILE *fe, FILE *svg,double *param,int *pext,int smooth)
{
  int fa=0,fs=1;
  double ang_start,ang_end;

  if (smooth) 
    {
      param[0] *= 1.25; param[1] *= 1.25;
      param[2] *= 1.25; param[3] *= 1.25;
      pext[4] *= 1.25; pext[5] *= 1.25;
      pext[6] *= 1.25; pext[7] *= 1.25;
    }
  double a = pow((param[2]-param[3])/(param[2]+param[3]),2);
  double C = M_PI*(param[2]+param[3])*(1+3*a/(10+sqrt(4-3*a)));

  ang_start = atan2(pext[5]-param[1],pext[4]-param[0]); 
  ang_end = atan2(pext[7]-param[1],pext[6]-param[0]); 
  if (angle_diff(ang_start,ang_end)<M_2__PI*SQRT2/C && angle_diff_signed(ang_start,ang_end)>0) 
    {
      fprintf(svg,"<ellipse transform=\"translate(%f %f) rotate(%f)\" rx=\"%f\" ry=\"%f\" stroke=\"red\" fill=\"none\" stroke-width= \"%d\" /> \n",param[0],param[1],(param[4])*180/M_PI,param[2],param[3],1);    
    }
  else 
    {
      double x1,y1,x2,y2;
      if (ang_start<0) ang_start+=M_2__PI;
      if (ang_end<0) ang_end+=M_2__PI;
  
      if (ang_end<ang_start) ang_end += M_2__PI;

      if ((ang_end-ang_start)>M_PI) fa = 1;

      rosin_point(param,pext[4],pext[5],&x1,&y1);

      rosin_point(param,pext[6],pext[7],&x2,&y2);

     
      fprintf(svg,"<path d=\"M %f,%f A%f,%f %f %d,%d %f,%f\"",
	x1,y1,param[2],param[3],param[4]*180/M_PI,fa,fs,x2,y2);
      fprintf(svg," fill=\"none\" stroke =\"red\" stroke-width=\"%d\" />\n",1);   
      //fprintf(svg," fill=\"none\" stroke =\"black\" stroke-width=\"%d\" />\n",1); 
      
    }
  fprintf(fe,"%f %f %f %f %f %f %f\n",param[0],param[1],param[2],param[3],param[4],ang_start,ang_end);

}
/*---------------------------------------------------------------------------*/
