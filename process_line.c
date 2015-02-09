/*----------------------------------------------------------------------------

  ELSD - Ellipse and Line Segment Detector

  Copyright (c) 2007-2011 rafael grompone von gioi (grompone@gmail.com)
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
#include <float.h>
#include "elsd.h"
#include "process_curve.h"



/*----------------------------------------------------------------------------*/
/** Copy one rectangle structure to another.
 */
void rect_copy(struct rect * in, struct rect * out)
{
  if( in == NULL || out == NULL ) error("rect_copy: invalid 'in' or 'out'.");
  out->x1 = in->x1;
  out->y1 = in->y1;
  out->x2 = in->x2;
  out->y2 = in->y2;
  out->width = in->width;
  out->x = in->x;
  out->y = in->y;
  out->theta = in->theta;
  out->dx = in->dx;
  out->dy = in->dy;
  out->prec = in->prec;
  out->p = in->p;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Rectangle points iterator auxiliary function.
 */
double inter_low(double x, double x1, double y1, double x2, double y2)
{
  if( x1 > x2 || x < x1 || x > x2 )
    {
      fprintf(stderr,"inter_low: x %g x1 %g x2 %g.\n",x,x1,x2);
      error("impossible situation.");
    }
  if( double_equal(x1,x2) && y1<y2 ) return y1;
  if( double_equal(x1,x2) && y1>y2 ) return y2;
  return y1 + (x-x1) * (y2-y1) / (x2-x1);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Rectangle points iterator auxiliary function.
 */
double inter_hi(double x, double x1, double y1, double x2, double y2)
{
  if( x1 > x2 || x < x1 || x > x2 )
    {
      fprintf(stderr,"inter_hi: x %g x1 %g x2 %g.\n",x,x1,x2);
      error("impossible situation.");
    }
  if( double_equal(x1,x2) && y1<y2 ) return y2;
  if( double_equal(x1,x2) && y1>y2 ) return y1;
  return y1 + (x-x1) * (y2-y1) / (x2-x1);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Free memory used by a rectangle iterator.
 */
void ri_del(rect_iter * iter)
{
  if( iter == NULL ) error("ri_del: NULL iterator.");
  free( (void *) iter );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Check if the iterator finished the full iteration.
 */
int ri_end(rect_iter * i)
{
  if( i == NULL ) error("ri_end: NULL iterator.");
  return (double)(i->x) > i->vx[2];
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Increment a rectangle iterator.
 */
void ri_inc(rect_iter * i)
{
  if( i == NULL ) error("ri_inc: NULL iterator.");

  if( (double) (i->x) <= i->vx[2] ) i->y++;

  while( (double) (i->y) > i->ye && (double) (i->x) <= i->vx[2] )
    {
      /* new x */
      i->x++;

      if( (double) (i->x) > i->vx[2] ) return; /* end of iteration */

      /* update lower y limit for the line */
      if( (double) i->x < i->vx[3] )
        i->ys = inter_low((double)i->x,i->vx[0],i->vy[0],i->vx[3],i->vy[3]);
      else i->ys = inter_low((double)i->x,i->vx[3],i->vy[3],i->vx[2],i->vy[2]);

      /* update upper y limit for the line */
      if( (double)i->x < i->vx[1] )
        i->ye = inter_hi((double)i->x,i->vx[0],i->vy[0],i->vx[1],i->vy[1]);
      else i->ye = inter_hi((double)i->x,i->vx[1],i->vy[1],i->vx[2],i->vy[2]);

      /* new y */
      i->y = (int) ceil(i->ys);
    }
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create and initialize a rectangle iterator.
 */
rect_iter * ri_ini(struct rect * r)
{
  double vx[4],vy[4];
  int n,offset;
  rect_iter * i;

  if( r == NULL ) error("ri_ini: invalid rectangle.");

  i = (rect_iter *) malloc(sizeof(rect_iter));
  if( i == NULL ) error("ri_ini: Not enough memory.");

  vx[0] = r->x1 - r->dy * r->width / 2.0;
  vy[0] = r->y1 + r->dx * r->width / 2.0;
  vx[1] = r->x2 - r->dy * r->width / 2.0;
  vy[1] = r->y2 + r->dx * r->width / 2.0;
  vx[2] = r->x2 + r->dy * r->width / 2.0;
  vy[2] = r->y2 - r->dx * r->width / 2.0;
  vx[3] = r->x1 + r->dy * r->width / 2.0;
  vy[3] = r->y1 - r->dx * r->width / 2.0;

  if( r->x1 < r->x2 && r->y1 <= r->y2 ) offset = 0;
  else if( r->x1 >= r->x2 && r->y1 < r->y2 ) offset = 1;
  else if( r->x1 > r->x2 && r->y1 >= r->y2 ) offset = 2;
  else offset = 3;

  for(n=0; n<4; n++)
    {
      i->vx[n] = vx[(offset+n)%4];
      i->vy[n] = vy[(offset+n)%4];
    }

  /* starting point */
  i->x = (int) ceil(i->vx[0]) - 1;
  i->y = (int) ceil(i->vy[0]);
  i->ys = i->ye = -DBL_MAX;

  /* advance to the first point */
  ri_inc(i);

  return i;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Compute a rectangle's NFA value.
 */
double rect_nfa(struct rect *rec,image_double angles,double logNT)
{
  rect_iter * i;
  int pts = 0;
  int alg = 0;
  if( rec == NULL ) error("rect_nfa: invalid rectangle.");
  if( angles == NULL ) error("rect_nfa: invalid 'angles'.");
  double theta = rec->theta - M_PI/2.0;
  for(i=ri_ini(rec); !ri_end(i); ri_inc(i))
    if( i->x >= 0 && i->y >= 0 &&
        i->x < (int) angles->xsize && i->y < (int) angles->ysize )
      {
        ++pts;
        if( isaligned(i->x, i->y, angles, theta, rec->prec) ) 
          {
	    ++alg;
	  }

      }
  ri_del(i);
  return nfa(pts,alg,rec->p,logNT);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*---------------------------------- Regions ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Compute a region's angle.
 */
double get_theta(struct point *reg,int reg_size0,int reg_size,double x,double y,
                         image_double modgrad,double reg_angle,double prec)
{
  double lambda1,lambda2,tmp,theta,weight,sum;
  double Ixx = 0.0;
  double Iyy = 0.0;
  double Ixy = 0.0;
  int i;
  

  /* check parameters */
  if( reg == NULL ) error("get_theta: invalid region.");
  if( reg_size <= 1 ) error("get_theta: region size <= 1.");
  if( modgrad == NULL || modgrad->data == NULL )
    error("get_theta: invalid 'modgrad'.");
  if( prec < 0.0 ) error("get_theta: 'prec' must be positive.");

  /*----------- theta ---------------------------------------------------*/
  /*
      Region inertia matrix A:
         Ixx Ixy
         Ixy Iyy
      where
        Ixx = \sum_i y_i^2
        Iyy = \sum_i x_i^2
        Ixy = -\sum_i x_i y_i

      lambda1 and lambda2 are the eigenvalues, with lambda1 >= lambda2.
      They are found by solving the characteristic polynomial
      det(\lambda I - A) = 0.

      To get the line segment direction we want to get the eigenvector of
      the smaller eigenvalue. We have to solve a,b in:
        a.Ixx + b.Ixy = a.lambda2
        a.Ixy + b.Iyy = b.lambda2
      We want the angle theta = atan(b/a). I can be computed with
      any of the two equations:
        theta = atan( (lambda2-Ixx) / Ixy )
      or
        theta = atan( Ixy / (lambda2-Iyy) )

      When |Ixx| > |Iyy| we use the first, otherwise the second
      (just to get better numeric precision).
   */
  sum = 0.0;
  for(i=reg_size0; i<reg_size; i++)
    {
      weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
      Ixx += ( (double) reg[i].y - y ) * ( (double) reg[i].y - y ) * weight;
      Iyy += ( (double) reg[i].x - x ) * ( (double) reg[i].x - x ) * weight;
      Ixy -= ( (double) reg[i].x - x ) * ( (double) reg[i].y - y ) * weight;
      sum += weight;
    }
  if( sum <= 0.0 ) error("get_theta: weights sum less or equal to zero.");
  Ixx /= sum;
  Iyy /= sum;
  Ixy /= sum;
  lambda1 = ( Ixx + Iyy + sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) ) / 2.0;
  lambda2 = ( Ixx + Iyy - sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) ) / 2.0;
  if( fabs(lambda1) < fabs(lambda2) )
    {
      fprintf(stderr,"Ixx %g Iyy %g Ixy %g lamb1 %g lamb2 %g - lamb1 < lamb2\n",
                      Ixx,Iyy,Ixy,lambda1,lambda2);
      tmp = lambda1;
      lambda1 = lambda2;
      lambda2 = tmp;
    }

  if( fabs(Ixx) > fabs(Iyy) )
    theta = atan2( lambda2-Ixx, Ixy );
  else
    theta = atan2( Ixy, lambda2-Iyy );

  /* The previous procedure doesn't care about orientation,
     so it could be wrong by 180 degrees. Here is corrected if necessary. */
  if( angle_diff(theta-M_1_2_PI,reg_angle) > prec ) theta += M_PI;

  if( angle_diff(theta-M_1_2_PI,reg_angle) > prec ) theta = reg_angle + M_1_2_PI;

  while (theta>M_PI) theta -= M_2__PI;
  while (theta<-M_PI) theta += M_2__PI;

  return theta;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Computes a rectangle that covers a region of points.
 */
void region2rect(struct point *reg,int reg_size0,int reg_size,image_double modgrad,
                 double reg_angle, double prec, double p, struct rect * rec)
{
  double x,y,dx,dy,l,w,theta,weight,sum,l_min,l_max,w_min,w_max;
  int i;

  /* check parameters */
  if( reg == NULL ) error("region2rect: invalid region.");
  if( reg_size <= 1 ) error("region2rect: region size <= 1.");
  if( modgrad == NULL || modgrad->data == NULL )
    error("region2rect: invalid image 'modgrad'.");
  if( rec == NULL ) error("region2rect: invalid 'rec'.");

  /* center */
  x = y = sum = 0.0;
  for(i=reg_size0; i<reg_size; i++)
    {
      weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
      x += (double) reg[i].x * weight;
      y += (double) reg[i].y * weight;
      sum += weight;
    }
  if( sum <= 0.0 ) error("region2rect: weights sum equal to zero.");
  x /= sum;
  y /= sum;

  /* theta */
  
  theta = get_theta(reg,reg_size0,reg_size,x,y,modgrad,reg_angle,prec);

  /* length and width */
  dx = cos(theta);
  dy = sin(theta);
  l_min = l_max = w_min = w_max = 0.0;
  for(i=reg_size0; i<reg_size; i++)
    {
      l =  ( (double) reg[i].x - x) * dx + ( (double) reg[i].y - y) * dy;
      w = -( (double) reg[i].x - x) * dy + ( (double) reg[i].y - y) * dx;

      if( l > l_max ) l_max = l;
      if( l < l_min ) l_min = l;
      if( w > w_max ) w_max = w;
      if( w < w_min ) w_min = w;
    }

  /* store values */
  rec->x1 = x + l_min * dx;
  rec->y1 = y + l_min * dy;
  rec->x2 = x + l_max * dx;
  rec->y2 = y + l_max * dy;
  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<1)
    {
      rec->x1 += (rec->dx+rec->dy); rec->y1 += (rec->dx+rec->dy);
      rec->x2 -= (rec->dx+rec->dy); rec->y2 -= (rec->dx+rec->dy);
    }
  rec->width = w_max - w_min;
  rec->x = x;
  rec->y = y;
  rec->theta = theta;
  rec->dx = dx;
  rec->dy = dy;
  rec->prec = prec;
  rec->p = p;

  if( rec->width < 1.0 ) 
    rec->width = 1.0;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Try some rectangles variations to improve NFA value.
    Only if the rectangle is not meaningful (i.e., log_nfa <= eps).
 */
static double rect_improve( struct rect * rec, image_double angles,
                            double logNT, double eps )
{
  struct rect r;
  double log_nfa,log_nfa_new;
  double delta = 0.5;
  double delta_2 = delta / 2.0;
  int n;

  log_nfa = rect_nfa(rec,angles,logNT);
 
  if( log_nfa > eps ) return log_nfa;
  
  /* try finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
    {
      r.p /= 2.0;
      r.prec = r.p * M_PI;
      log_nfa_new = rect_nfa(&r,angles,logNT);
      if( log_nfa_new > log_nfa )
        {
          log_nfa = log_nfa_new;
          rect_copy(&r,rec);
        }
    }

  if( log_nfa > eps ) return log_nfa;

  /* try to reduce width */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
    {
      if( (r.width - delta) >= 0.5 )
        {
          r.width -= delta;
          log_nfa_new = rect_nfa(&r,angles,logNT);
          if( log_nfa_new > log_nfa )
            {
              rect_copy(&r,rec);
              log_nfa = log_nfa_new;
            }
        }
    }

  if( log_nfa > eps ) return log_nfa;

  /* try to reduce one side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
    {
      if( (r.width - delta) >= 0.5 )
        {
          r.x1 += -r.dy * delta_2;
          r.y1 +=  r.dx * delta_2;
          r.x2 += -r.dy * delta_2;
          r.y2 +=  r.dx * delta_2;
          r.width -= delta;
          log_nfa_new = rect_nfa(&r,angles,logNT);
          if( log_nfa_new > log_nfa )
            {
              rect_copy(&r,rec);
              log_nfa = log_nfa_new;
            }
        }
    }

  if( log_nfa > eps ) return log_nfa;

  /* try to reduce the other side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
    {
      if( (r.width - delta) >= 0.5 )
        {
          r.x1 -= -r.dy * delta_2;
          r.y1 -=  r.dx * delta_2;
          r.x2 -= -r.dy * delta_2;
          r.y2 -=  r.dx * delta_2;
          r.width -= delta;
          log_nfa_new = rect_nfa(&r,angles,logNT);
          if( log_nfa_new > log_nfa )
            {
              rect_copy(&r,rec);
              log_nfa = log_nfa_new;
            }
        }
    }

  if( log_nfa > eps ) return log_nfa;

  /* try even finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
    {
      r.p /= 2.0;
      r.prec = r.p * M_PI;
      log_nfa_new = rect_nfa(&r,angles,logNT);
      if( log_nfa_new > log_nfa )
        {
          log_nfa = log_nfa_new;
          rect_copy(&r,rec);
        }
    }

  return log_nfa;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Reduce the region size, by elimination the points far from the
    starting point, until that leads to rectangle with the right
    density of region points or to discard the region if too small.
 */
int reduce_region_radius(struct point * reg, int reg_size0, int * reg_size,
                         image_double modgrad, image_double gradx,
                         image_double grady, double reg_angle,
                         double prec, double p,struct rect * rec,
                         image_char used, image_double angles,
                         double density_th, double logNT, double eps)
{
  double density,rad1,rad2,rad,xc,yc,log_nfa;
  int i;
  /* check parameters */
  if( reg == NULL ) error("refine: invalid pointer 'reg'.");
  if( reg_size == NULL ) error("refine: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("refine: 'prec' must be positive.");
  if( rec == NULL ) error("refine: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("refine: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("refine: invalid image 'angles'.");

  /* compute region points density */
  density = (double) (*reg_size - reg_size0) /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
  //printf("red rad dens init %f \n",density);
  if( density >= density_th ) return TRUE;

  /* compute region radius */
  xc = (double) reg[reg_size0].x;
  yc = (double) reg[reg_size0].y;
  rad1 = dist( xc, yc, rec->x1, rec->y1 );
  rad2 = dist( xc, yc, rec->x2, rec->y2 );
  rad = rad1 > rad2 ? rad1 : rad2;

  while( density < density_th )
    {
      rad *= 0.75;

      /* remove points from the region and update 'used' map */
      for(i=reg_size0; i<*reg_size; i++)
        if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) > rad )
          {
            /* point not kept, mark it as NOTUSED */
            used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
            /* remove point from the region */
            reg[i].x = reg[*reg_size-1].x; /* if i==*reg_size-1 copy itself */
            reg[i].y = reg[*reg_size-1].y;
            --(*reg_size);
            --i; /* to avoid skipping one point */
          }
      /* reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */

      if( (*reg_size - reg_size0) < 2 ) return FALSE;

      /* re-compute rectangle */
      region2rect(reg,reg_size0,*reg_size,modgrad,reg_angle,prec,p,rec);
      if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<1)
        {
          rec->x1 += (rec->dx+rec->dy); rec->y1 += (rec->dx+rec->dy);
          rec->x2 -= (rec->dx+rec->dy); rec->y2 -= (rec->dx+rec->dy);
        }

      /* try to improve the rectangle and compute NFA */
      log_nfa = rect_improve(rec,angles,logNT,eps);
      /* re-compute region points density */
      density = (double) (*reg_size - reg_size0) /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );     
    }

  /* if the final rectangle is meaningful accept, otherwise reject */
  if( log_nfa > eps ) return TRUE;
  else return FALSE;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Refine a rectangle. For that, an estimation of the angle tolerance is
    performed by the standard deviation of the angle at points near the
    region's starting point. Then, a new region is grown starting from the
    same point, but using the estimated angle tolerance.
    If this fails to produce a rectangle with the right density of
    region points, 'reduce_region_radius' is called to try to
    satisfy this condition.
 */
int refine(struct point *reg,int reg_size0,int *reg_size,image_double modgrad,
           image_double gradx,image_double grady,
           double prec,double p,struct rect *rec,
           image_char used,image_double angles,double density_th,
           double logNT,double eps)
{
  double angle,ang_d,mean_angle,tau,density,xc,yc,ang_c,sum,s_sum,log_nfa;
  double reg_angle;
  int i,n;
 
  /* check parameters */
  if( reg == NULL ) error("refine: invalid pointer 'reg'.");
  if( reg_size == NULL ) error("refine: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("refine: 'prec' must be positive.");
  if( rec == NULL ) error("refine: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("refine: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("refine: invalid image 'angles'.");

  /* compute region points density */
  density = (double) (*reg_size-reg_size0) /
            ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
 
  if( density >= density_th ) return TRUE;

  /*------ First try: reduce angle tolerance ------*/

  /* compute the new mean angle and tolerance */
  xc = (double) reg[reg_size0].x;
  yc = (double) reg[reg_size0].y;
  ang_c = angles->data[ reg[reg_size0].x + reg[reg_size0].y * angles->xsize ];
  sum = s_sum = 0.0;
  n = 0;
  for(i=reg_size0; i<*reg_size; i++)
    {
      used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
      if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) < rec->width )
        {
          angle = angles->data[ reg[i].x + reg[i].y * angles->xsize ];
          ang_d = angle_diff_signed(angle,ang_c);
          sum += ang_d;
          s_sum += ang_d * ang_d;
          ++n;
        }
    }
  mean_angle = sum / (double) n;
  tau = max(2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                         + mean_angle*mean_angle ),0.1); /* 2 * standard deviation */
  /*tau = 2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                         + mean_angle*mean_angle );*/ /* 2 * standard deviation */

 
  /* find a new region from the same starting point and new angle tolerance */
  *reg_size = reg_size0 + 1;
  region_grow(reg[reg_size0].x,reg[reg_size0].y,angles,reg,reg_size,&reg_angle,used,tau);

  /* if the region is too small, reject */
  if( *reg_size-reg_size0 <= 2 ) return FALSE;

  /* re-compute rectangle */
  region2rect(reg,reg_size0,*reg_size,modgrad,reg_angle,prec,p,rec);
  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<1)
    {
      rec->x1 += (rec->dx+rec->dy); rec->y1 += (rec->dx+rec->dy);
      rec->x2 -= (rec->dx+rec->dy); rec->y2 -= (rec->dx+rec->dy);
    }

  /* try to improve the rectangle and compute NFA */
  log_nfa = rect_improve(rec,angles,logNT,eps);

  /* re-compute region points density */
  density = (double) (*reg_size - reg_size0) /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radius(reg,reg_size0,reg_size,modgrad,gradx,grady,reg_angle,prec,p,
                                rec,used,angles,density_th,logNT,eps);

  /* if the final rectangle is meaningful accept, otherwise reject */
  if( log_nfa > eps ) return TRUE;
  else return FALSE;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Decide if a given rectangle is a valid line segment.  
 */
void valid_line(struct point* reg,int *reg_size,double reg_angle,double prec,
                   double p,struct rect *rec,double *l,image_double grad,
                   image_double gradx,image_double grady, 
                   image_char used,image_double angles,double density_th, 
                   double logNT,double mlog10eps,double *nfaline)
{
  *nfaline = mlog10eps;
  l[4] = 0;
  if (*reg_size>2)
    {      
      region2rect(reg,0,*reg_size,grad,reg_angle,prec,p,rec);
      if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<1)
        {
          rec->x1 += (rec->dx+rec->dy); rec->y1 += (rec->dx+rec->dy);
          rec->x2 -= (rec->dx+rec->dy); rec->y2 -= (rec->dx+rec->dy);
        }
      if(!refine(reg,0,reg_size,grad,gradx,grady,prec,p,rec,used,angles,
            density_th, logNT, mlog10eps)) 
        {
          l[4]=0;
        }
      else
        { 
          *nfaline = rect_improve(rec,angles,logNT,mlog10eps);
          if((*nfaline) <= mlog10eps) l[4]=0;
          else
            {
              l[0] = rec->x1 + 0.5;
              l[1] = rec->y1 + 0.5;
              l[2] = rec->x2 + 0.5;
              l[3] = rec->y2 + 0.5;
              l[4] = rec->width;
            }
          rec->x1 += 0.5; rec->y1 += 0.5;
          rec->x2 += 0.5; rec->y2 += 0.5;
        }
    }
}
/*----------------------------------------------------------------------------*/

