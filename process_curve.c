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
//#include <f2c.h>
//#include <clapack.h>



/*---------------------------------------------------------------------------*/
/** Convert ellipse from matrix form to common form:
    ellipse = (centrex,centrey,ax,ay,orientation).
 */
void ellipse2param(double *p,double param[])
{
	/* p = [ a,     1/2*b, 1/2*d; 
	         1/2*b,   c,   1/2*e;
	         1/2*d, 1/2*e,   f   ]; */
  double a,b,c,d,e,f;
  a = p[0];
  b = 2*p[1];
  c = p[4];
  d = 2*p[2];
  e = 2*p[5];
  f = p[8]; 
  double thetarad,cost,sint,cos_squared,sin_squared,cos_sin,Ao,Au,Av,Auu,Avv,
  tuCentre,tvCentre,wCentre,uCentre,vCentre,Ru,Rv;

  thetarad=0.5*atan2(b,a-c); 
  cost=cos(thetarad);
  sint=sin(thetarad);
  sin_squared=sint*sint;
  cos_squared=cost*cost;
  cos_sin=sint*cost;
  Ao=f;
  Au=d*cost+e* sint;
  Av=-d*sint+e* cost;
  Auu=a*cos_squared+c*sin_squared+b*cos_sin;
  Avv=a*sin_squared+c*cos_squared-b*cos_sin;

  if(Auu==0 || Avv==0){ param[0]=0;param[1]=0;param[2]=0;param[3]=0;param[4]=0;}
  else
    {
      tuCentre=-Au/(2.*Auu);
      tvCentre=-Av/(2.*Avv);
      wCentre=Ao-Auu*tuCentre*tuCentre-Avv*tvCentre*tvCentre;
      uCentre=tuCentre*cost-tvCentre*sint;
      vCentre=tuCentre*sint+tvCentre*cost;
      Ru=-wCentre/Auu;
      Rv=-wCentre/Avv;
      if (Ru>0) Ru=pow(Ru,0.5);
      else Ru=-pow(-Ru,0.5);
      if (Rv>0) Rv=pow(Rv,0.5);
      else Rv=-pow(-Rv,0.5);
      param[0]=uCentre;param[1]=vCentre;
      param[2]=Ru;param[3]=Rv;param[4]=thetarad;
    }
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Compute matrix mean. 
    'sz' x 'dim' = matrix size (lines x cols)
    'm' = output mean (vector of size 1 x 'dim').
 */ 
void mean(double *v,int sz,double *m,int dim)
{
  int i,j;
  double sum[dim];
  if (v == NULL) error("mean: Invalid pointer");
  if (sz<=0) error("mean: Invalid size");
  if (dim<=0) error("mean: Invalid size");
  for (i=0;i<dim;i++) sum[i] = 0.0;
  for (i=0;i<sz;i++) 
    for (j=0;j<dim;j++)    
      sum[j]+=v[i*dim+j];

  for (i=0;i<dim;i++) m[i] = sum[i]/(double)sz;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Compute normalisation matrix (translates and normalises a set of 
    2D homogeneous points so that their centroid is at the origin and 
    their mean distance from the origin is sqrt(2).
    Input points in 'v' are in cartesian coordinates.
    Output matrix 'T' is 3 x 3.
 */
void vgg_conditioner_from_points(double *T, struct point *reg, int sz)
{
  /* check parameters */
  if (reg == NULL) error("vgg_conditioner_from_points: invalid points list");
  if (sz<=0) error("vgg_conditioner_from_points: invalid size list");
  double m[2] = {0.0, 0.0};
  double Qt = 0.0;
  double Qmean = 0.0, val = 0.0;
  int i;
  
  /* compute mean point */
  for (i=0; i<sz; i++)
    {
      m[0] += reg[i].x;
      m[1] += reg[i].y;
    }  
  m[0] /= (double)sz; m[1] /= (double)sz;

  /* compute variance */
  for (i=0; i<sz; i++)
    {
      val = (reg[i].x - m[0])*(reg[i].x - m[0])+(reg[i].y - m[1])*(reg[i].y - m[1]);
      Qt += sqrt(val);
    }  	
  Qmean = Qt/(double)sz; 

  val = SQRT2/Qmean;

  T[1] = T[3] = 0;
  T[0] = T[4] = val;
  T[2] = -val*m[0];
  T[5] = -val*m[1];
  T[6] = T[7] = 0;
  T[8] = 1;    
}  
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** antisym(u)  A = [ 0,-u(3),u(2); u(3),0,-u(1); -u(2),u(1),0 ]; 
 */
void antisym(double *u,double *A)
{
  A[0] = A[4] = A[8] = 0;
  A[1] = -u[2];
  A[2] = u[1];
  A[3] = u[2];
  A[5] = -u[0];
  A[6] = -u[1];
  A[7] = u[0];
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Compute equations for circle/ellipse fitting.
 */
void fit_equations(struct point *reg,int reg_size,image_double gradx, 
                   image_double grady, double *vgg)
{
  /* check parameters */
  if(reg == NULL) error("fit_equations: invalid region");
  if (reg_size<=0) error("fit_equations: invalid region size");
  if (reg_size*24>gSizeBufferDouble)
    {
      gBufferDouble = (double*)realloc(gBufferDouble, sizeof(double) * reg_size * 24);  
      if (!gBufferDouble) error("fit_equations: not enough memory");
      gSizeBufferDouble = reg_size*24;
    }
  int i,j;
  double K[27];
  double asym[9];
  int idx;
  double crosspr[3]; 
  double pnormx, pnormy, dirnormx, dirnormy;
  int addr;

  /* compute normalisation matrix */
  vgg_conditioner_from_points(vgg,reg,reg_size);

  /* compute equation system */
  for (i=0;i<reg_size;i++)
    { 
      idx = i*4*6;
      /* normalise point (pnormx,pnormy) = VGG*(x,y) */
      pnormx = vgg[0]*reg[i].x+vgg[1]*reg[i].y+vgg[2];
      pnormy = vgg[3]*reg[i].x+vgg[4]*reg[i].y+vgg[5];
      
      /* normalise gradient direction (dirnormx,dirnormy) = VGG*(dx,dy) */
      addr = reg[i].y*gradx->xsize+reg[i].x;
      dirnormx = -vgg[0] * grady->data[addr] + vgg[1] * gradx->data[addr];
      dirnormy = -vgg[3] * grady->data[addr] + vgg[4] * gradx->data[addr];

      /* cross product (pnormx,pnormy) x (dirnormx,dirnormy) = tangent line */
      crosspr[0] = - dirnormy;
      crosspr[1] = dirnormx;
      crosspr[2] = pnormx * dirnormy - pnormy * dirnormx;

      /* tangent's equation : eq = -transpose(kron(TPts(1:3,i),antisym(l)))*J; */
      antisym(crosspr,asym);

      for (j=0;j<9;j++) K[j] = asym[j]*pnormx; 
      for (j=0;j<9;j++) K[j+9] = asym[j]*pnormy; 
      for (j=0;j<9;j++) K[j+18] = asym[j];

      gBufferDouble[idx]   = -K[0];            gBufferDouble[idx+6]   = -K[1];            gBufferDouble[idx+12]   = -K[2];
      gBufferDouble[idx+1] = -(K[3]+K[9]);     gBufferDouble[idx+6+1] = -(K[4]+K[10]);    gBufferDouble[idx+12+1] = -(K[5]+K[11]); 
      gBufferDouble[idx+2] = -(K[6]+K[18]);    gBufferDouble[idx+6+2] = -(K[7]+K[19]);    gBufferDouble[idx+12+2] = -(K[8]+K[20]);
      gBufferDouble[idx+3] = -K[12];           gBufferDouble[idx+6+3] = -K[13];           gBufferDouble[idx+12+3] = -K[14]; 
      gBufferDouble[idx+4] = -(K[15]+K[21]);   gBufferDouble[idx+6+4] = -(K[16]+K[22]);   gBufferDouble[idx+12+4] = -(K[17]+K[23]);
      gBufferDouble[idx+5] = -K[24];           gBufferDouble[idx+6+5] = -K[25];           gBufferDouble[idx+12+5] = -K[26]; 

      /* position's equation : eq = transpose(kron(TPts(1:3,i),TPts(1:3,i)))*J; */
      gBufferDouble[idx+18] = pnormx * pnormx; 
      gBufferDouble[idx+19] = 2 * pnormx * pnormy; 
      gBufferDouble[idx+20] = 2 * pnormx;
      gBufferDouble[idx+21] = pnormy * pnormy; 
      gBufferDouble[idx+22] = 2 * pnormy; 
      gBufferDouble[idx+23] = 1;
    }
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Algebraic circle fitting using positional and tangential constraints.
 */
void fitcircle(int reg_size, double *vgg, double *param)
{
  /* check parameters */
  if (reg_size<=0) error("fitcircle: invalid size");
  int i,j,k;
  double A[16]; 
  int idx;
  
  idx = 0;

  for(i=0; i<reg_size*4*6; i+=6)
    { 
      gBufferDouble[i] = gBufferDouble[i] + gBufferDouble[i+3];
      gBufferDouble[i+1] = gBufferDouble[i+2];
      gBufferDouble[i+2] = gBufferDouble[i+4];
      gBufferDouble[i+3] = gBufferDouble[i+5];
    } 
  /* A = EQ'*EQ; */
  for (i=0;i<16;i++) A[i] = 0.0;

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4*reg_size;k++)
        A[i*4+j] += gBufferDouble[k*6+i]*gBufferDouble[k*6+j];

#define SIZE4 4
  char JOBZ = 'V';
  char UPLO = 'U';
  int M = SIZE4;
  int LDA = M;
  int LWORK = 4*SIZE4;
  int INFO;
  double W[SIZE4];
  double WORK[LWORK];
  dsyev_(&JOBZ, &UPLO, &M, A, &LDA, W, WORK, &LWORK, &INFO); 
 
  double s[9];
  s[0] = s[4] = A[0];
  s[1] = s[3] = 0;
  s[2] = s[6] = A[1];
  s[5] = s[7] = A[2];
  s[8] = A[3];

  /* apply inverse(normalisation matrix) */  
  /* C = T'*[ x(1),0,x(2); 0,x(1),x(3) ; x(2),x(3),x(4)]*T;*/

  double C[9];
  C[0] = vgg[0]*vgg[0]*s[0]+vgg[3]*vgg[3]*s[4];
  C[1] = vgg[0]*vgg[1]*s[0]+vgg[3]*vgg[4]*s[4];
  C[2] = vgg[0]*vgg[2]*s[0]+vgg[3]*vgg[5]*s[4]+vgg[0]*s[2]+vgg[3]*s[5];

  C[3] = vgg[0]*vgg[1]*s[0]+vgg[3]*vgg[4]*s[4];
  C[4] = vgg[1]*vgg[1]*s[0]+vgg[4]*vgg[4]*s[4];
  C[5] = vgg[1]*vgg[2]*s[0]+vgg[4]*vgg[5]*s[4]+vgg[1]*s[2]+vgg[4]*s[5];

  C[6] = vgg[0]*vgg[2]*s[0]+vgg[0]*s[6]+vgg[3]*vgg[5]*s[4]+vgg[3]*s[7];
  C[7] = vgg[1]*vgg[2]*s[0]+vgg[1]*s[6]+vgg[4]*vgg[5]*s[4]+vgg[4]*s[7];
  C[8] = vgg[2]*vgg[2]*s[0]+vgg[2]*s[6]+vgg[5]*vgg[5]*s[4]+vgg[5]*s[7]+vgg[2]*s[2]+vgg[5]*s[5]+s[8];
  
  ellipse2param(C, param);
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Algebraic ellipse fitting using positional and tangential constraints.
 */
void fitellipse(int reg_size, double *vgg, double *param)
{
  /* check parameters */
  if (reg_size<=0) error("fitellipse: invalid size");
  int i,j,k;
  double A[36];   
 
  /* A = EQ'*EQ; */
  for (i=0;i<36;i++) A[i] = 0.0;

  for (i=0;i<6;i++)
    for (j=0;j<6;j++)
      for (k=0;k<4*reg_size;k++)
        A[i*6+j] += gBufferDouble[k*6+i]*gBufferDouble[k*6+j];

#define SIZE6 6
  char JOBZ = 'V';
  char UPLO = 'U';    
  int M = SIZE6;
  int LDA = M;
  int LWORK = 4*SIZE6;
  int INFO;
  double W[SIZE6];
  double WORK[LWORK];
  dsyev_(&JOBZ, &UPLO, &M, A, &LDA, W, WORK, &LWORK, &INFO); 

    
  double s[9];
  s[0] = A[0];
  s[1] = s[3] = A[1];
  s[2] = s[6] = A[2];
  s[4] = A[3];
  s[5] = s[7] = A[4];
  s[8] = A[5];
  
  /* apply inverse(normalisation matrix) */
  /* C = T'*[ x(1),x(2),x(3); x(2),x(4),x(5) ; x(3),x(5),x(6)]*T; */
  
  double C[9];
  C[0] = vgg[0]*vgg[0]*s[0]+vgg[0]*vgg[3]*s[3]+vgg[0]*vgg[3]*s[1]+vgg[3]*vgg[3]*s[4];
  C[1] = vgg[0]*vgg[1]*s[0]+vgg[1]*vgg[3]*s[3]+vgg[0]*vgg[4]*s[1]+vgg[3]*vgg[4]*s[4];
  C[2] = vgg[0]*vgg[2]*s[0]+vgg[2]*vgg[3]*s[3]+vgg[0]*vgg[5]*s[1]+vgg[3]*vgg[5]*s[4]+vgg[0]*s[2]+vgg[3]*s[5];

  C[3] = vgg[0]*vgg[1]*s[0]+vgg[0]*vgg[4]*s[3]+vgg[1]*vgg[3]*s[1]+vgg[3]*vgg[4]*s[4];
  C[4] = vgg[1]*vgg[1]*s[0]+vgg[1]*vgg[4]*s[3]+vgg[1]*vgg[4]*s[1]+vgg[4]*vgg[4]*s[4];
  C[5] = vgg[1]*vgg[2]*s[0]+vgg[2]*vgg[4]*s[3]+vgg[1]*vgg[5]*s[1]+vgg[4]*vgg[5]*s[4]+vgg[1]*s[2]+vgg[4]*s[5];

  C[6] = vgg[0]*vgg[2]*s[0]+vgg[0]*vgg[5]*s[3]+vgg[0]*s[6]+vgg[2]*vgg[3]*s[1]+vgg[3]*vgg[5]*s[4]+vgg[3]*s[7];
  C[7] = vgg[1]*vgg[2]*s[0]+vgg[1]*vgg[5]*s[3]+vgg[1]*s[6]+vgg[2]*vgg[4]*s[1]+vgg[4]*vgg[5]*s[4]+vgg[4]*s[7];
  C[8] = vgg[2]*vgg[2]*s[0]+vgg[2]*vgg[5]*s[3]+vgg[2]*s[6]+vgg[2]*vgg[5]*s[1]+vgg[5]*vgg[5]*s[4]+vgg[5]*s[7]+vgg[2]*s[2]+vgg[5]*s[5]+s[8];
  
  ellipse2param(C, param);
}
/*---------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute the min value of a double array.
 */
double min_array(double *v, int sz)
{
  double m = DBL_MAX;
  int i;
  for (i=0;i<sz;i++) 
    if (v[i]<m) m = v[i];
  return m;
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** Compute the position of the min value in a double array.
 */
void min_array_pos(double *v, int sz,int *pos)
{
  double m = DBL_MAX;
  int i;
  for (i=0;i<sz;i++) 
    if (v[i]<m) 
      {
        m = v[i];
        *pos = i;
      }
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute the point belonging to an ellipse, closest to a given point with
    integer coordinates.
 */
void rosin_point (double *param,double x,double y,double *xi,double *yi)
{  

  double ae2 = param[2]*param[2];
  double be2 = param[3]*param[3];
  x = x - param[0];
  y = y - param[1];
  double xp = x*cos(-param[4])-y*sin(-param[4]);
  double yp = x*sin(-param[4])+y*cos(-param[4]);
  double fe2;
  fe2 = ae2-be2;
  double X = xp*xp;
  double Y = yp*yp;
  double delta = (X+Y+fe2)*(X+Y+fe2)-4*X*fe2;
  double A = (X + Y + fe2 - sqrt(delta))/2.0; 
  double ah = sqrt(A);
  double bh2 = fe2-A;
  double term = (A*be2+ae2*bh2);
  double xx = ah*sqrt(ae2*(be2+bh2)/term);
  double yy = param[3]*sqrt(bh2*(ae2-A)/term);
  double d[4];
  int pos;

  d[0] = dist(xp,yp,xx,yy);
  d[1] = dist(xp,yp,xx,-yy);
  d[2] = dist(xp,yp,-xx,yy);
  d[3] = dist(xp,yp,-xx,-yy);
  min_array_pos(d,4,&pos);
  switch (pos)
    { 
      case 0: break;
      case 1: yy = -yy; break;
      case 2: xx = -xx; break;
      case 3: xx = -xx; yy = -yy; break;
      default: break;
    }
  *xi = xx*cos(param[4])-yy*sin(param[4])+param[0];
  *yi = xx*sin(param[4])+yy*cos(param[4])+param[1];
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Approximate the distance between a point and an ellipse using Rosin distance.
 */
double d_rosin (double *param, double x, double y)
{ 
  double ae2 = param[2]*param[2];
  double be2 = param[3]*param[3];
  x = x - param[0];
  y = y - param[1];
  double xp = x*cos(-param[4])-y*sin(-param[4]);
  double yp = x*sin(-param[4])+y*cos(-param[4]);
  double fe2;
  fe2 = ae2-be2;
  double X = xp*xp;
  double Y = yp*yp;
  double delta = (X+Y+fe2)*(X+Y+fe2)-4*X*fe2;
  double A = (X + Y + fe2 - sqrt(delta))/2.0; 
  double ah = sqrt(A);
  double bh2 = fe2-A;
  double term = (A*be2+ae2*bh2);
  double xi = ah*sqrt(ae2*(be2+bh2)/term);
  double yi = param[3]*sqrt(bh2*(ae2-A)/term);
  double d[4];


  d[0] = dist(xp,yp,xi,yi);
  d[1] = dist(xp,yp,xi,-yi);
  d[2] = dist(xp,yp,-xi,yi);
  d[3] = dist(xp,yp,-xi,-yi);
  double dmin = min_array(d,4); 
  
  if (X+Y>xi*xi+yi*yi)
    return dmin;
  else return -dmin; 
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute the angle of a point belonging to an ellipse using the focal property
 */
double angle(double x, double y, double *foci)
{
  double tmp1 = atan2(y-foci[1], x-foci[0]);
  double tmp2 = atan2(y-foci[3], x-foci[2]);
  double theta;
  double tmp3 = angle_diff_signed(tmp1,tmp2);
  theta = tmp1-tmp3/2.0;
  while( theta <= -M_PI ) theta += M_2__PI;
  while( theta >   M_PI ) theta -= M_2__PI;
  return theta;
}
/*----------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
int in_interval(double x,double a,double b){
if (x>=a&&x<b) return 1;
return 0;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Determine if the gradient converges or diverges to/from the centre.
 */
int int_ext(double xc, double yc, int px, int py, image_double angles)
{
  double a;
  int dir = 0;
  a = atan2(py-yc,px-xc);
  if (angle_diff(a,angles->data[py*angles->xsize+px])<M_1_2_PI)
    dir = 1;
  return dir;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Return +1 for positive value, and -1 for strictly negative value.
 */
int sign(double val)
{
  if (val>=0) return 1;
  else return -1;
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/** Compute the pixel seed for new region growing.
 */
void px_seed(struct rect rec, struct point *reg, int reg_size, int reg_size0, 
             image_char used, image_double angles, image_double grad, int *pext, 
             int *ps, double prec)
{
  int xx,yy,i;
  int adr;
  int xtmp = pext[0], ytmp = pext[1];
  double lmax = 0.0, l;
  double gradmax = 0;

  /* find the most extreme point of a region */
  for (i = reg_size0;i<reg_size;i++)
    {
      l = fabs(((double)reg[i].x-xtmp)*rec.dx +((double)reg[i].y-ytmp)*rec.dy);
      if (l>lmax)
        {
          pext[0] = reg[i].x;
          pext[1] = reg[i].y;
          lmax = l;
        }
    }

  /* find a pixel with strong gradient, near the end of the current region */
  ps[0] = pext[0]; ps[1] = pext[1];

  for (xx=pext[0]-1;xx<=pext[0]+1;xx++)
    for(yy = pext[1]-1;yy<=pext[1]+1;yy++)
      {
        if (xx>=0 && xx<(int)used->xsize && yy >=0 && yy< (int)used->ysize)
          {
            adr = yy*used->xsize+xx;
            if (used->data[adr] != USED && angle_diff(angles->data[adr],rec.theta-M_1_2_PI)<2*prec)
              if (grad->data[adr]>gradmax && dist(xx,yy,xtmp,ytmp)>=lmax)
                {
                  ps[0] = xx;
                  ps[1] = yy;
                  gradmax = grad->data[adr];
                }
          }
      }    
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/** Compute polygonal approximation of a curve.
 */
void subcurve(struct point *reg, int *reg_size, double prec, double p,
              image_double angles, image_char used, image_double grad,
              image_double gradx, image_double grady, double logNT, 
              double mlog10eps, double density_th, int *pext, int *ps, 
              int *sgn, double *ang0, double *spir)
{
  int convex = 1;
  int reg_size0 = 0;
  struct rect rec;
  double reg_angle;
  int flag;
  int i;
  double difang;

  init_rect(&rec);
  while (*reg_size != reg_size0 && convex)
    {
      reg_size0 = *reg_size;
      ++(*reg_size);
      region_grow(ps[0], ps[1], angles, reg, reg_size, &reg_angle, used, prec);   
          
      if (*reg_size-reg_size0 > 1)
        {
          region2rect(reg, reg_size0, *reg_size, grad, reg_angle, prec, p, &rec);
          flag = refine(reg, reg_size0, reg_size, grad, gradx, grady, prec, p, &rec, used,
                        angles, density_th, logNT, mlog10eps);             
          if (dist(rec.x1,rec.y1,rec.x2,rec.y2)<1)
            {
              rec.x1 += (rec.dx+rec.dy); rec.y1 += (rec.dx+rec.dy);
              rec.x2 -= (rec.dx+rec.dy); rec.y2 -= (rec.dx+rec.dy);
            }
          
          if (*sgn == 0) 
            *sgn = sign(angle_diff_signed(rec.theta, *ang0));  
          
          if (*reg_size-reg_size0 > 1)
            {
              difang = angle_diff_signed(rec.theta, *ang0);
              /* convexity & smoothness check & non-spiral check */
              if (sign(difang) == *sgn && fabs(difang) < M_1_2_PI && 
                 (*spir = *spir + fabs(difang)) <= M_2__PI)
                { 
                  *ang0 = rec.theta;
                  px_seed(rec, reg, *reg_size, reg_size0, used, angles, grad, pext, ps, prec);   
                }      
              else /* not convex or not smooth */
                {
                  convex = 0;
                  /* clean last rectangle */
                  for (i=reg_size0;i<*reg_size;i++) 
                    used->data[reg[i].y*used->xsize+reg[i].x] = NOTUSED;   
                  /* substract the contribution of the last rectangle */
                  *reg_size = reg_size0;
                }
            }               
          else 
            { 
              convex = 0;
              *reg_size = reg_size0;
            }
        }
      else 
        { 
          convex = 0;
          *reg_size = reg_size0;
        }      
    }
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
void curve_grow(struct point *reg, int *reg_size, struct rect rec, double prec, 
                double p, image_double angles, image_char used, image_double grad,
                image_double gradx, image_double grady, double logNT, 
                double mlog10eps, double density_th)
{
  int ps1[2], ps2[2];
  int pext1[2], pext2[2];
  int sgn = 0;
  double ang0 = rec.theta;
  double spir = 0.0;

  /* compute right and left pixel seed using the ends of the initial rectangle */
  pext2[0] = rec.x1; pext2[1] = rec.y1;
  px_seed(rec, reg, *reg_size, 0, used, angles, grad, pext2, ps2, prec); 

  pext1[0] = rec.x2; pext1[1] = rec.y2; 
  px_seed(rec, reg, *reg_size, 0, used, angles, grad, pext1, ps1, prec);
  
  /* scan in one direction */ 
  subcurve(reg, reg_size, prec, p, angles, used, grad, gradx, grady, logNT, 
           mlog10eps, density_th, pext2, ps2, &sgn, &ang0, &spir);

  /* scan the other direction if a complete tour was not scanned */ 
  ang0 = rec.theta;
  sgn = -sgn;
  if (spir < M_2__PI)
    subcurve(reg, reg_size, prec, p, angles, used, grad, gradx, grady, logNT, 
             mlog10eps, density_th, pext1, ps1, &sgn, &ang0, &spir);
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
void valid_curve(struct point *reg,int *reg_size,double prec,double p,
                 image_double angles,image_char used,image_double grad,
                 image_double gradx,image_double grady,double *paramc,
                 double *parame,struct rect rec,double *logNT,double mlog10eps,
                 double density_th,int *min_size,double *nfa,int *pext,
                 struct point3 *regc, struct point3 *rege, int *regp_size)
{
  /* check params */
  if (reg == NULL) error("valid_curve : null region");
  if (angles == NULL || angles->data == NULL) error("valid_curve : null angles");
  if (grad == NULL || grad->data == NULL) error("valid_curve : null grad");
  if (gradx == NULL || gradx->data == NULL) error("valid_curve : null gradx");
  if (grady == NULL || grady->data == NULL) error("valid_curve : null grady");
  if (used == NULL || used->data == NULL) error("valid_curve : null used");
  if (*reg_size <= 0) error("valid_curve : invalid reg_size");
  
  int ell_ok;
  double vgg[9];
  int dir;
  nfa[1] = nfa[2] = mlog10eps; 
  regp_size[1] = regp_size[2] = 0;
  curve_grow(reg,reg_size,rec,prec,p,angles,used,grad,gradx,grady,logNT[0],mlog10eps,density_th);
  if (*reg_size>min_size[1])
    {
      fit_equations(reg,*reg_size,gradx,grady,vgg);
      /* perform first ellipse fitting, and then circle fitting because circle fitting modifies the equations 'eq' */
      if (*reg_size>min_size[2])
        {
          fitellipse(*reg_size,vgg,parame);
        }
      fitcircle(*reg_size,vgg,paramc);  
      dir = int_ext(paramc[0], paramc[1], reg[0].x, reg[0].y, angles);
      
      nfa[1] = valid_circle(reg,*reg_size,used,prec,p,angles,grad,gradx,grady,paramc, 
                            logNT[1],dir,pext,regc,&regp_size[1],min_size[1],mlog10eps);

      if (*reg_size>min_size[2])
        {
          ell_ok = check_ellipse(parame);

          if (ell_ok)
            nfa[2] = valid_ellipse(reg,*reg_size,used,prec,p,angles,grad,gradx,grady,parame, 
                                   logNT[2],dir,pext,rege,&regp_size[2],min_size[2],mlog10eps);
        }
    }
}
/*---------------------------------------------------------------------------*/

