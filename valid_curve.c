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
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include "elsd.h"
#include "valid_curve.h"
#include "process_curve.h"
#include "process_line.h"


/*---------------------------------------------------------------------------*/
/** Compute max element in an array and return the max value and its position.
 */
double max_array(double *a, int sz, int *poz)
{
  /* check parameters */
  if (a == NULL) error("max: invalid pointer");
  int i;
  double max = (double)LONG_MIN;
  for (i=0;i<sz;i++)
    if (max<a[i]) {max = a[i]; *poz = i;}
  return max;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Quicksort the elements of an array and return the sorted array 
    and the original positions of the elements.
 */
void quickSort(double *arr, int elements, int *pos) {

  #define  MAX_LEVELS  300

  int  beg[MAX_LEVELS], end[MAX_LEVELS], i, L, R, swap ;
  double piv;
  for (i=0;i<elements;i++) pos[i] = i;
  int ptmp;
  beg[0]=0; end[0]=elements;
  i = 0;
  while (i>=0) 
    {
      L=beg[i]; R=end[i]-1;
      if (L<R) 
        {
          piv=arr[L];
          ptmp = pos[L];
          while (L<R) 
            {
              while (arr[R]>=piv && L<R) R--; 
              if (L<R) 
                {
                  arr[L]=arr[R]; pos[L] = pos[R]; L++;
                }
              while (arr[L]<=piv && L<R) L++; 
              if (L<R) 
                {
                  arr[R]=arr[L]; pos[R] = pos[L]; R--;
                }
            }
          arr[L]=piv; pos[L]=ptmp; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;

          if (end[i]-beg[i]>end[i-1]-beg[i-1]) 
            {
              swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
              swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; 
            }
        }
      else
        {
          i--; 
        }
    }
}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute the delimiting angles of a circular/elliptical arc. Return 
    the delimiting angles and the index positions of the extreme points.
 */
void extreme_sorted_angles(int sz, double *ang_start,double *ang_end, int *idx)
{ 
  /* check parameters */
  if (sz<=0) error("extreme_angles : invalid size of angles list");
  int i;
  double alphamax = 0.0;
  int poz;
  double difftmp;

  if (sz>gSizeBufferInt)
    {
      gBufferInt = (int*)realloc(gBufferInt, sz*sizeof(int));  
      if (!gBufferInt) error("extreme_sorted_angles: not enough memory");
      gSizeBufferInt = sz;
    }

  quickSort(gBufferDouble,sz,gBufferInt);

  for (i=0;i<sz-1;i++)
    {
      difftmp = gBufferDouble[i+1]-gBufferDouble[i];
      if (difftmp>alphamax)
        {
          alphamax = difftmp;
          poz = i;
        }
    }
  difftmp = gBufferDouble[0] + M_2__PI - gBufferDouble[sz-1];
  if (difftmp>alphamax)
    {
      alphamax = difftmp;
      poz = sz - 1;
    }
  *ang_end = gBufferDouble[poz];
  idx[1] = gBufferInt[poz];
  if (poz != (sz-1)) 
    {
      *ang_start = gBufferDouble[poz+1]; 
      idx[0] = gBufferInt[poz+1];
    }
  else 
    {
      *ang_start = gBufferDouble[0]; 
      idx[0] = gBufferInt[0];
    }
}
/*---------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
/** Test if angle 'ang' is in interval ['ang1', 'ang2']; 
    'ang', 'ang1', 'ang2' are in [0 , 2pi].
 */
int isInAng(double ang, double ang1, double ang2)
{
  
  int ok = 0;
  if (ang2>ang1) 
    {
      if (ang >= ang1 && ang <= ang2) ok = 1;
    }
  else
    if (ang >= ang1 || ang <= ang2) ok = 1;
  return ok;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Clean region.
 */
void clean_reg(struct point* reg, int reg_size, image_char used)
{
  int i = 0;
  for (i=0;i<reg_size;i++) 
    used->data[reg[i].y*used->xsize+reg[i].x] = NOTUSED;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute a circle's NFA value : discrete NFA.
 */
double valid_circle(struct point *reg, int reg_size, image_char used, double prec, 
                    double p, image_double angles, image_double grad, image_double gradx,
                    image_double grady, double *param, double logNTC, int dir, int *pext,
                    struct point3 *regc, int *regc_size, int min_size, double mlog10eps)
{
  int alg = 0;
  double d_min = angles->xsize;
  double d_max = -d_min;
  double d,theta,theta0;
  int i,xx,yy;  
  double ang_start,ang_end;
  int idx[2];
  
  /* check parameters */ 
  if( angles == NULL ) error("valid_circle: invalid 'angles'.");

  /* clean region for validation */
  clean_reg(reg,reg_size,used);
  

  /* refine: search for connected aligned points on the given circle starting from the same seed point */
  used->data[reg[0].y*used->xsize+reg[0].x] = USED;
  *regc_size = 1;
  for(i=0; i<*regc_size; i++)
    for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
      for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
        if (xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
            used->data[xx+yy*used->xsize] != USED)
          { 
            theta0 = atan2((double)yy-param[1],(double)xx-param[0]);

            if (dir==0)
              {
	        if (theta0>0) theta = -(M_PI-theta0);
		else theta = M_PI + theta0;			
	      }
            else theta = theta0;
            if (isaligned( xx, yy, angles, theta, prec)) 
              {  
                used->data[xx+yy*used->xsize] = USED;
                reg[*regc_size].x = xx;
                reg[*regc_size].y = yy;                    
                ++(*regc_size);
              }  
          }

  if (*regc_size>min_size)
    {
      /* reestimate circle on the connected aligned points to have a better precision */
      double vgg[9];
      fit_equations(reg,*regc_size,gradx,grady,vgg);      
      fitcircle(*regc_size,vgg,param);

      /* compute parameters of the circular ring: width and delimiting angles */
      for(i=0;i<*regc_size;i++)
        {
          /* compute width and store angles in global temp gArray1 */
          d = sqrt((reg[i].x-param[0])*(reg[i].x-param[0])+(reg[i].y-param[1])*(reg[i].y-param[1]))-param[2];
          if (d<d_min) d_min = d;
          if (d>d_max) d_max = d;
          gBufferDouble[i] = atan2((double)reg[i].y-param[1],(double)reg[i].x-param[0]);
          if (gBufferDouble[i]<0) gBufferDouble[i] += M_2__PI;
          used->data[reg[i].y*grad->xsize+reg[i].x] = NOTUSED;
        }

      /* compute delimiting angles */
      if (*regc_size>2) extreme_sorted_angles(*regc_size,&ang_start,&ang_end,idx); 

      /* extract extreme contour points */
      pext[0] = reg[idx[0]].x; pext[1] = reg[idx[0]].y; 
      pext[2] = reg[idx[1]].x; pext[3] = reg[idx[1]].y;

      /* scan the circular ring and count the number of points and the number of aligned points */
      *regc_size = 1;
      used->data[reg[0].y*used->xsize+reg[0].x] = USEDCIRC;
      regc[0].x = reg[0].x;
      regc[0].y = reg[0].y;
      regc[0].z = USEDCIRC;
      for(i=0; i<*regc_size; i++)
        for(xx=regc[i].x-1; xx<=regc[i].x+1; xx++)
          for(yy=regc[i].y-1; yy<=regc[i].y+1; yy++)
            if (xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
                used->data[xx+yy*used->xsize] != USED && used->data[xx+yy*used->xsize] != USEDCIRC)
              { 
                d = sqrt((xx-param[0])*(xx-param[0])+(yy-param[1])*(yy-param[1]))-param[2];
                theta0 = atan2((double)yy-param[1],(double)xx-param[0]);
                if (theta0<0) theta = theta0 + M_2__PI;
                else theta = theta0;

                if(d>=d_min && d<=d_max && isInAng(theta,ang_start,ang_end))
                  {
                    if (dir==0)
                      {
		        if (theta0>0) theta = -(M_PI-theta0);
		        else theta = M_PI + theta0;			
	              }
                    else theta = theta0;
                    if (isaligned( xx, yy, angles, theta, prec)) 
                      {  
                        ++alg;
                        regc[*regc_size].z = USEDCIRC;
                      }
                    else 
                      {
                        regc[*regc_size].z = USEDCIRCNA;
                      }
                    used->data[xx+yy*used->xsize] = USEDCIRC;
                    regc[*regc_size].x = xx;
                    regc[*regc_size].y = yy;
                    
                    ++(*regc_size);
                  }
              }
      return nfa(*regc_size,alg,p,logNTC); /* compute NFA value */ 
    }
  else
    return mlog10eps;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute ellipse foci, given ellipse params.
 */
void ellipse_foci(double *param, double *foci)
{
  double f = sqrt(param[2]*param[2]-param[3]*param[3]);
  foci[0] = param[0]+f*cos(param[4]);
  foci[1] = param[1]+f*sin(param[4]);
  foci[2] = param[0]-f*cos(param[4]);
  foci[3] = param[1]-f*sin(param[4]);
}
/*----------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Compute an ellipse's NFA value : discrete NFA.
 */
double valid_ellipse(struct point *reg, int reg_size, image_char used, double prec,
                     double p, image_double angles, image_double grad, image_double gradx,
                     image_double grady, double *param, double logNTE, int dir, int *pext, 
                     struct point3 *rege, int *rege_size, int min_size, double mlog10eps)
{
  int alg = 0;
  double d_min = angles->xsize;
  double d_max = -d_min;
  double d,theta,theta0;
  int i,xx,yy; 
  double ang_start,ang_end;
  int idx[2];
  double foci[4]; /* xf1, yf1, xf2, yf2 */
 
  /* check parameters */ 
  if( angles == NULL ) error("valid_ellipse: invalid 'angles'.");

  ellipse_foci(param, foci);

  /* refine: search for connected aligned points on the given ellipse starting from the same seed point */
  *rege_size = 1;
  used->data[reg[0].y*used->xsize+reg[0].x] = USED;
  for(i=0; i<*rege_size; i++)
    for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
      for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
        if (xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
            used->data[xx+yy*used->xsize] != USED)
          { 
            theta = angle((double)xx, (double)yy, foci);
            if (dir==0)
              {
		if (theta>0) theta = -(M_PI-theta);
		else theta = M_PI + theta;			
	      }
            if (isaligned( xx, yy, angles, theta, prec)) 
              {
                used->data[xx+yy*used->xsize] = USED;      
                reg[*rege_size].x = xx;
                reg[*rege_size].y = yy;
                ++(*rege_size);
              }
          }

  if (*rege_size>min_size)
    {
      /* reestimate ellipse on the connected aligned points to have a better precision */
      double vgg[9];
      fit_equations(reg,*rege_size,gradx,grady,vgg);
      fitellipse(*rege_size,vgg,param);
      int ell_ok = check_ellipse(param);
      if (ell_ok)
        {
          ellipse_foci(param, foci);

          /* compute parameters of the elliptical ring: width and delimiting angles */
          for(i=0;i<*rege_size;i++)
            {
              /* compute width and store angles in global temp gArray1 */
              d = d_rosin(param,(double)reg[i].x, (double)reg[i].y);
              if (d<d_min) d_min = d;
              if (d>d_max) d_max = d;
              gBufferDouble[i] = atan2((double)reg[i].y-param[1],(double)reg[i].x-param[0]);
              if (gBufferDouble[i]<0) gBufferDouble[i] += M_2__PI;
              used->data[reg[i].y*grad->xsize+reg[i].x] = NOTUSED;
            }

          /* compute delimiting angles */
          if (*rege_size>2) extreme_sorted_angles(*rege_size, &ang_start, &ang_end, idx); 

          /* extract extreme contour points */
          pext[4] = reg[idx[0]].x; pext[5] = reg[idx[0]].y; 
          pext[6] = reg[idx[1]].x; pext[7] = reg[idx[1]].y;

          /* scan the elliptical ring and count the number of points and the number of aligned points */
          *rege_size = 1;
          used->data[reg[0].y*used->xsize+reg[0].x] = USEDELL;
          rege[0].x = reg[0].x;
          rege[0].y = reg[0].y;
          rege[0].z = USEDELL;
          
          for(i=0; i<*rege_size; i++)
            for(xx=rege[i].x-1; xx<=rege[i].x+1; xx++)
              for(yy=rege[i].y-1; yy<=rege[i].y+1; yy++)
                if (xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
                used->data[xx+yy*used->xsize] != USED && used->data[xx+yy*used->xsize] != USEDELL)
                  { 
                    d = d_rosin(param,(double)xx,(double)yy);
                    theta0 = atan2((double)yy-param[1],(double)xx-param[0]);
                    if (theta0<0) theta = theta0 + M_2__PI;
                    else theta = theta0;
                    if(d>=d_min && d<=d_max && isInAng(theta,ang_start,ang_end))
                      {
                        theta = angle((double)xx, (double)yy, foci);
                        if (dir==0)
                          {
		            if (theta>0) theta = -(M_PI-theta);
		            else theta = M_PI + theta;			
		          }
                        if (isaligned( xx, yy, angles, theta, prec)) 
                          {
                            ++alg;
                            rege[*rege_size].z = USEDELL;
                          }
                        else 
                          { 
                            rege[*rege_size].z = USEDELLNA;
                          }
                        used->data[xx+yy*used->xsize] = USEDELL;
                        rege[*rege_size].x = xx;
                        rege[*rege_size].y = yy;
                        ++(*rege_size);
                      }
                  }
          return nfa(*rege_size,alg,p,logNTE); /* compute NFA value */ 
        }
      else 
        { 
          *rege_size = 0;
          return mlog10eps;
        }
    }
  else
    return mlog10eps;
}
/*----------------------------------------------------------------------------*/
