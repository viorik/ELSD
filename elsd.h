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

#ifndef ELSD_H
#define ELSD_H


#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

#define NOTDEF -1024.0
#define M_3_2_PI 4.71238898038
#define M_1_2_PI 1.57079632679
#define M_2__PI  6.28318530718
#define NOTUSED 0
#define USED 1
#define USEDCIRC 2
#define USEDCIRCNA 3
#define USEDELL 4
#define USEDELLNA 5

#define SQRT2 1.414213562373095

#define max(A, B) (((A)>(B))?(A):(B))
#define min(A, B) (((A)<(B))?(A):(B))



/*----------------------------------------------------------------------------*/
/* -------------------------- Global temporary variables -------------------- */
/*----------------------------------------------------------------------------*/
extern double *gBufferDouble;
extern int *gBufferInt;
extern int gSizeBufferDouble,gSizeBufferInt;
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*--------------------------- Rectangle structure ----------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
struct rect /* line segment with width */
{
  double x1,y1,x2,y2;  /* first and second point of the line segment */
  double width;        /* rectangle width */
  double x,y;          /* center of the rectangle */
  double theta;        /* angle */
  double dx,dy;        /* vector with the line segment angle */
  double prec;         /* tolerance angle */
  double p;            /* probability of a point with angle within 'prec' */
};
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
typedef struct rit 
{
  double vx[4];
  double vy[4];
  double ys,ye;
  int x,y;
} rect_iter;
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
typedef struct gauss_filter_s
{
int dim;
double sigma;
double mean;
double *values;
} *gauss_filter;
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
typedef struct image_char_s
{
  unsigned char *data;
  unsigned int xsize,ysize;
} *image_char;

void free_image_char(image_char i);
image_char new_image_char(unsigned int xsize, unsigned int ysize);
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value );
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** int image data type
 */

typedef struct image_int_s
{
  int * data;
  unsigned int xsize,ysize;
} *image_int;

void free_image_int(image_int i);
image_int new_image_int(unsigned int xsize, unsigned int ysize);
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value );
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** double image data type
 */
typedef struct image_double_s
{
  double *data;
  unsigned int xsize,ysize;
} *image_double;

void free_image_double(image_double i);
image_double new_image_double(unsigned int xsize, unsigned int ysize);
image_double new_image_double_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value );
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
struct coorlist
{
  int x,y;
  struct coorlist * next;
};
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
struct point {int x,y;};
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
struct point3 {int x,y,z;};
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
struct doublepoint {double x,y;};
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void error(char * msg);
int double_equal(double a, double b);
int isaligned(int x, int y, image_double theta, double angle,
                     double precision);
double dist(double x1, double y1, double x2, double y2);

void get_cadran(double angle, int* cad);

int check_ellipse(double *param);
double nfa(int n, int k, double p, double logNT);
void region_grow(int x, int y, image_double angles, struct point * reg,
                 int * reg_size, double * reg_angle, image_char used,
                 double prec );
double angle_diff_signed(double a, double b);
double angle_diff(double a, double b);
double angle_diff_full(double a, double b, int sens);
/*----------------------------------------------------------------------------*/
#endif

