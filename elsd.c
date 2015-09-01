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
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "elsd.h"
#include "write_svg.h"
#include "valid_curve.h"
#include "process_curve.h"
#include "process_line.h"

/*----------------------------------------------------------------------------*/
/* Init global temporary variables */
int gSizeBufferDouble = 1;
int gSizeBufferInt    = 1;

double *gBufferDouble; 
int    *gBufferInt;    
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Euclidean distance between two points (x1,y1) and (x2,y2).
 */
double dist(double x1, double y1, double x2, double y2)
{
    return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}
/*----------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
void error(char * msg)
{
    fprintf(stderr,"%s\n",msg);
    exit(EXIT_FAILURE);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
double angle_diff_full(double a, double b, int sens)
{
  a -= b;
  if (sens == 1 && a<0) a+= M_2__PI;
  if (sens == 2)
    if (a>0) 
      a = M_2__PI - a;
    else 
      a = -a;
  return a;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
double angle_diff(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  if( a < 0.0 ) a = -a;
  return a;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Signed angle difference.
 */
double angle_diff_signed(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  return a;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Check two double numbers for equality. 
 */
#define RELATIVE_ERROR_FACTOR 100.0
int double_equal(double a, double b)
{
    double abs_diff,aa,bb,abs_max;

    if( a == b ) return TRUE;

    abs_diff = fabs(a-b);
    aa = fabs(a);
    bb = fabs(b);
    abs_max = aa > bb ? aa : bb;

    /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
    if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

    /* equal if relative error <= factor x eps */
    return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*------------------------------ PGM image I/O -------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments(FILE * f)
{
  int c;
  do
    {
      while(isspace(c=getc(f))); /* skip spaces */
      if(c=='#') /* skip comments */
        while( c!='\n' && c!='\r' && c!=EOF )
          c=getc(f);
    }
  while( c == '#' || isspace(c) );
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Read an ASCII number from a PGM file.
 */
static unsigned int get_num(FILE * f)
{
  unsigned int num;
  int c;

  while(isspace(c=getc(f)));
  if(!isdigit(c)) error("Error: corrupted PGM file.");
  num = (unsigned int) (c - '0');
  while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");

  return num;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_double".
    If the name is "-" the file is read from standard input.
 */
static image_double read_pgm_image_double(char * name)
{
  FILE * f;
  int c,bin;
  unsigned int xsize,ysize,depth,x,y;
  image_double image;

  /* open file */
  f = fopen(name,"rb");
  if( f == NULL ) error("Error: unable to open input image file.");

  /* read header */
  if( getc(f) != 'P' ) error("Error: not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("Error: not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);            /* X size */
  skip_whites_and_comments(f);
  ysize = get_num(f);            /* Y size */
  skip_whites_and_comments(f);
  depth = get_num(f);            /* depth */
  if(depth==0) fprintf(stderr,"Warning: depth=0, probably invalid PGM file\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("Error: corrupted PGM file.");

  /* get memory */
  image = new_image_double(xsize,ysize);

  /* read data */
  for(y=0;y<ysize;y++)
    for(x=0;x<xsize;x++)
      image->data[ x + y * xsize ] = bin ? (double) getc(f)
                                         : (double) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
    error("Error: unable to close file while reading PGM file.");

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*----------------------------- Image Data Types -----------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Free memory used in image_char 'i'.
 */
void free_image_char(image_char i)
{
  if( i == NULL || i->data == NULL )
    error("free_image_char: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize'.
 */
image_char new_image_char(unsigned int xsize, unsigned int ysize)
{
  image_char image;

  /* check parameters */
  if( xsize == 0 || ysize == 0 ) error("new_image_char: invalid image size.");

  /* get memory */
  image = (image_char) malloc( sizeof(struct image_char_s) );
  if( image == NULL ) error("not enough memory.");
  image->data = (unsigned char *) calloc( (size_t) (xsize*ysize),
                                          sizeof(unsigned char) );
  if( image->data == NULL ) error("not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value )
{
  image_char image = new_image_char(xsize,ysize); /* create image */
  unsigned int N = xsize*ysize;
  unsigned int i;

  /* check parameters */
  if( image == NULL || image->data == NULL )
    error("new_image_char_ini: invalid image.");

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = fill_value;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Free memory used in image_int 'i'.
 */
void free_image_int(image_int i)
{
  if( i == NULL || i->data == NULL )
    error("free_image_int: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize'.
 */
image_int new_image_int(unsigned int xsize, unsigned int ysize)
{
  image_int image;

  /* check parameters */
  if( xsize == 0 || ysize == 0 ) error("new_image_int: invalid image size.");

  /* get memory */
  image = (image_int) malloc( sizeof(struct image_int_s) );
  if( image == NULL ) error("not enough memory.");
  image->data = (int *) calloc( (size_t) (xsize*ysize), sizeof(int) );
  if( image->data == NULL ) error("not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value )
{
  image_int image = new_image_int(xsize,ysize); /* create image */
  unsigned int N = xsize*ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = fill_value;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Free memory used in image_double 'i'.
 */
void free_image_double(image_double i)
{
  if( i == NULL || i->data == NULL )
    error("free_image_double: invalid input image.");
  free( (void *) i->data );
  free( (void *) i );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize'.
 */
image_double new_image_double(unsigned int xsize, unsigned int ysize)
{
  image_double image;

  /* check parameters */
  if( xsize == 0 || ysize == 0 ) error("new_image_double: invalid image size.");

  /* get memory */
  image = (image_double) malloc( sizeof(struct image_double_s) );
  if( image == NULL ) error("not enough memory.");
  image->data = (double *) calloc( (size_t) (xsize*ysize), sizeof(double) );
  if( image->data == NULL ) error("not enough memory.");

  /* set image size */
  image->xsize = xsize;
  image->ysize = ysize;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize',
    initialized to the value 'fill_value'.
 */
image_double new_image_double_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value )
{
  image_double image = new_image_double(xsize,ysize); /* create image */
  unsigned int N = xsize*ysize;
  unsigned int i;

  /* initialize */
  for(i=0; i<N; i++) image->data[i] = fill_value;

  return image;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*----------------------------- NFA computation ------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
   Calculates the natural logarithm of the absolute value of
   the gamma function of x using the Lanczos approximation,
   see http://www.rskey.org/gamma.htm

   The formula used is
     \Gamma(x) = \frac{ \sum_{n=0}^{N} q_n x^n }{ \Pi_{n=0}^{N} (x+n) }
                 (x+5.5)^(x+0.5) e^{-(x+5.5)}
   so
     \log\Gamma(x) = \log( \sum_{n=0}^{N} q_n x^n ) + (x+0.5) \log(x+5.5)
                     - (x+5.5) - \sum_{n=0}^{N} \log(x+n)
   and
     q0 = 75122.6331530
     q1 = 80916.6278952
     q2 = 36308.2951477
     q3 = 8687.24529705
     q4 = 1168.92649479
     q5 = 83.8676043424
     q6 = 2.50662827511
 */
double log_gamma_lanczos(double x)
{
  static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477,
                         8687.24529705, 1168.92649479, 83.8676043424,
                         2.50662827511 };
  double a = (x+0.5) * log(x+5.5) - (x+5.5);
  double b = 0.0;
  int n;

  for(n=0;n<7;n++)
    {
      a -= log( x + (double) n );
      b += q[n] * pow( x, (double) n );
    }
  return a + log(b);
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/**
   Calculates the natural logarithm of the absolute value of
   the gamma function of x using Robert H. Windschitl method,
   see http://www.rskey.org/gamma.htm

   The formula used is
     \Gamma(x) = \sqrt(\frac{2\pi}{x}) ( \frac{x}{e}
                   \sqrt{ x\sinh(1/x) + \frac{1}{810x^6} } )^x
   so
     \log\Gamma(x) = 0.5\log(2\pi) + (x-0.5)\log(x) - x
                     + 0.5x\log( x\sinh(1/x) + \frac{1}{810x^6} ).

   This formula is a good approximation when x > 15.
 */
double log_gamma_windschitl(double x)
{
  return 0.918938533204673 + (x-0.5)*log(x) - x
         + 0.5*x*log( x*sinh(1/x) + 1/(810.0*pow(x,6.0)) );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/**
   Calculates the natural logarithm of the absolute value of
   the gamma function of x. When x>15 use log_gamma_windschitl(),
   otherwise use log_gamma_lanczos().
 */
#define log_gamma(x) ((x)>15.0?log_gamma_windschitl(x):log_gamma_lanczos(x))

/*----------------------------------------------------------------------------*/
/*
   Computes -log10(NFA)

   NFA stands for Number of False Alarms:

     NFA = NT.b(n,k,p)

   NT    - number of tests
   b(,,) - tail of binomial distribution with parameters n,k and p

   The value -log10(NFA) is equivalent but more intuitive than NFA:
    -1 corresponds to 10 mean false alarms
     0 corresponds to 1 mean false alarm
     1 corresponds to 0.1 mean false alarms
     2 corresponds to 0.01 mean false alarms
     ...

   Used this way, the bigger the value, better the detection,
   and a logarithmic scale is used.

   Parameters:
     n,k,p - binomial parameters.
     logNT - logarithm of Number of Tests
 */
#define TABSIZE 100000
double nfa(int n, int k, double p, double logNT)
{
  static double inv[TABSIZE];   /* table to keep computed inverse values */
  double tolerance = 0.1;       /* an error of 10% in the result is accepted */
  double log1term,term,bin_term,mult_term,bin_tail,err,p_term;
  int i;
  //printf("%d %d \n", n,k);
  if( n<0 || k<0 || k>n || p<=0.0 || p>=1.0 )
    error("nfa: wrong n, k or p values.");

  if( n==0 || k==0 ) return -logNT;
  if( n==k ) return -logNT - (double) n * log10(p);

  p_term = p / (1.0-p);

  /* compute the first term of the series */
  /*
     binomial_tail(n,k,p) = sum_{i=k}^n bincoef(n,i) * p^i * (1-p)^{n-i}
     where bincoef(n,i) are the binomial coefficients.
     But
       bincoef(n,k) = gamma(n+1) / ( gamma(k+1) * gamma(n-k+1) ).
     We use this to compute the first term. Actually the log of it.
   */
  log1term = log_gamma( (double) n + 1.0 ) - log_gamma( (double) k + 1.0 )
           - log_gamma( (double) (n-k) + 1.0 )
           + (double) k * log(p) + (double) (n-k) * log(1.0-p);
  term = exp(log1term);

  /* in some cases no more computations are needed */
  if( double_equal(term,0.0) )        /* the first term is almost zero */
    {
      if( (double) k > (double) n * p )     /* at begin or end of the tail?  */
        return -log1term / M_LN10 - logNT;  /* end: use just the first term  */
      else
        return -logNT;                      /* begin: the tail is roughly 1  */
    }

  /* compute more terms if needed */
  bin_tail = term;
  for(i=k+1;i<=n;i++)
    {
      /*
         As
           term_i = bincoef(n,i) * p^i * (1-p)^(n-i)
         and
           bincoef(n,i)/bincoef(n,i-1) = n-1+1 / i,
         then,
           term_i / term_i-1 = (n-i+1)/i * p/(1-p)
         and
           term_i = term_i-1 * (n-i+1)/i * p/(1-p).
         1/i is stored in a table as they are computed,
         because divisions are expensive.
         p/(1-p) is computed only once and stored in 'p_term'.
       */
      bin_term = (double) (n-i+1) * ( i<TABSIZE ?
                   ( inv[i]!=0.0 ? inv[i] : ( inv[i] = 1.0 / (double) i ) ) :
                   1.0 / (double) i );

      mult_term = bin_term * p_term;
      term *= mult_term;
      bin_tail += term;
      if(bin_term<1.0)
        {
          /* When bin_term<1 then mult_term_j<mult_term_i for j>i.
             Then, the error on the binomial tail when truncated at
             the i term can be bounded by a geometric series of form
             term_i * sum mult_term_i^j.                            */
          err = term * ( ( 1.0 - pow( mult_term, (double) (n-i+1) ) ) /
                         (1.0-mult_term) - 1.0 );

          /* One wants an error at most of tolerance*final_result, or:
             tolerance * abs(-log10(bin_tail)-logNT).
             Now, the error that can be accepted on bin_tail is
             given by tolerance*final_result divided by the derivative
             of -log10(x) when x=bin_tail. that is:
             tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
             Finally, we truncate the tail if the error is less than:
             tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
          if( err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail ) break;
        }
    }
  //printf("nfa %d %d %f \n",n,k, -log10(bin_tail)- logNT);
  return -log10(bin_tail) - logNT;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Free memory used in a gauss_filter 'in'
 */
void free_gauss_filter(gauss_filter in)
{
  if( in == NULL || in->values == NULL )
    error("free_gauss_filter: invalid filter input.");
  free( (void *) in->values );
  free( (void *) in );
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Allocate space for a new gauss_filter of size 'dim'.
 */
gauss_filter new_gauss_filter(unsigned int dim)
{
  gauss_filter fil;

  if( dim <= 0 ) error("new_gauss_filter: invalid filter size");

  fil = (gauss_filter) malloc( sizeof(struct gauss_filter_s) );
  if( fil == NULL ) error("new_gauss_filter: not enough memory");
  fil->values = (double *) calloc( dim, sizeof(double) );
  if( fil->values == NULL ) error("new_gauss_filter: not enough memory");

  fil->dim = dim;
  fil->mean = 0.0;

  return fil;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/**
   Compute a Gaussian kernel of length 'kernel->dim',
   standard deviation 'sigma', and centered at value 'mean'.
   For example, if mean=0.5, the Gaussian will be centered
   in the middle point between values 'kernel->values[0]'
   and 'kernel->values[1]'.
 */
static void gaussian_kernel(gauss_filter kernel, double sigma, double mean)
{
  double sum = 0.0;
  double val;
  int i;

  if( kernel == NULL || kernel->values == NULL )
    error("gaussian_kernel: invalid struct 'kernel'.");
  if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");

  /* compute gaussian kernel */
  //kernel->size = 1;
  kernel->mean = mean;
  kernel->sigma = sigma; 
  for(i=0;i<kernel->dim;i++)
    {
      val = ( (double) i - mean ) / sigma;
      kernel->values[i] = exp( -0.5 * val * val );
      sum += kernel->values[i];
    }

  /* normalization */
  if( sum >= 0.0 ) for(i=0;i<kernel->dim;i++) kernel->values[i] /= sum;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/**
   Subsample image 'in' with Gaussian filtering, to a scale 'scale'
   (for example, 0.8 will give a result at 80% of the original size),
   using a standard deviation sigma given by:

     sigma = sigma_scale / scale,   if scale <  1.0
     sigma = sigma_scale,           if scale >= 1.0
 */
static image_double gaussian_sampler( image_double in, double scale,
                                      double sigma_scale )
{
  image_double aux;
  image_double out;
  gauss_filter kernel;
  unsigned int N,M,h,n,x,y;
  int i;
  int xc,yc,j,double_x_size,double_y_size;
  double sigma,xx,yy,sum,prec;
  
  if( in == NULL || in->data == NULL || in->xsize <= 0 || in->ysize <= 0 )
    error("gaussian_sampler: invalid image.");
  if( scale <= 0.0 ) error("gaussian_sampler: 'scale' must be positive.");
  if( sigma_scale <= 0.0 )
    error("gaussian_sampler: 'sigma_scale' must be positive.");

  /* get memory for images */
  N = (unsigned int) floor( in->xsize * scale );
  M = (unsigned int) floor( in->ysize * scale );
  aux = new_image_double(N,in->ysize);
  out = new_image_double(N,M);

  /* sigma, kernel size and memory for the kernel */
  sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
  /*
     The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  h = (unsigned int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1+2*h; /* kernel size */
  kernel = new_gauss_filter(n);

  /* auxiliary double image size variables */
  double_x_size = (int) (2 * in->xsize);
  double_y_size = (int) (2 * in->ysize);

  /* First subsampling: x axis */
  for(x=0;x<aux->xsize;x++)
    {
      /*
         x   is the coordinate in the new image.
         xx  is the corresponding x-value in the original size image.
         xc  is the integer value, the pixel coordinate of xx.
       */
      xx = (double) x / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with xc=0 get the values of xx from -0.5 to 0.5 */
      xc = (int) floor( xx + 0.5 );
      gaussian_kernel( kernel, sigma, (double) h + xx - (double) xc );
      /* the kernel must be computed for each x because the fine
         offset xx-xc is different in each case */

      for(y=0;y<aux->ysize;y++)
        {
          sum = 0.0;
          for(i=0;i<kernel->dim;i++)
            {
              j = xc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_x_size;
              while( j >= double_x_size ) j -= double_x_size;
              if( j >= (int) in->xsize ) j = double_x_size-1-j;

              sum += in->data[ j + y * in->xsize ] * kernel->values[i];
            }
          aux->data[ x + y * aux->xsize ] = sum;
        }
    }

  /* Second subsampling: y axis */
  for(y=0;y<out->ysize;y++)
    {
      /*
         y   is the coordinate in the new image.
         yy  is the corresponding x-value in the original size image.
         yc  is the integer value, the pixel coordinate of xx.
       */
      yy = (double) y / scale;
      /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with yc=0 get the values of yy from -0.5 to 0.5 */
      yc = (int) floor( yy + 0.5 );
      gaussian_kernel( kernel, sigma, (double) h + yy - (double) yc );
      /* the kernel must be computed for each y because the fine
         offset yy-yc is different in each case */

      for(x=0;x<out->xsize;x++)
        {
          sum = 0.0;
          for(i=0;i<kernel->dim;i++)
            {
              j = yc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_y_size;
              while( j >= double_y_size ) j -= double_y_size;
              if( j >= (int) in->ysize ) j = double_y_size-1-j;

              sum += aux->data[ x + j * aux->xsize ] * kernel->values[i];
            }
          out->data[ x + y * out->xsize ] = sum;
        }
    }

  /* free memory */
  free_gauss_filter(kernel);
  free_image_double(aux);
  return out;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Gradient orientation.
 */
image_double ll_angle(image_double in,double threshold,struct coorlist **list_p,
                      void **mem_p,image_double *gradx, image_double *grady, 
                      image_double *grad,unsigned int n_bins, double max_grad)
{
  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize <= 0 || in->ysize <= 0 )
    error("ll_angle: invalid image");
  if( threshold < 0.0 ) error("ll_angle: 'threshold' must be positive");
  if( list_p == NULL ) error("ll_angle: NULL pointer 'list_p'");
  if( mem_p == NULL ) error("ll_angle: NULL pointer 'mem_p'");
  if( n_bins <= 0 ) error("ll_angle: 'n_bins' must be positive");
  if( max_grad <= 0.0 ) error("ll_angle: 'max_grad' must be positive");
  image_double angles;
  unsigned int xsize,ysize,adr,ind,i,j;
  double com1,com2,gx,gy,norm,norm2;
  /* the rest of the variables are used for pseudo-ordering
     the gradient magnitude values */
  int list_count = 0;
  struct coorlist *list;
  struct coorlist **range_l_s; /* array of pointers to start of bin list */
  struct coorlist **range_l_e; /* array of pointers to end of bin list */
  struct coorlist *start;
  struct coorlist *end;


  xsize = in->xsize;
  ysize = in->ysize;
  /* allocate output image */
  angles = new_image_double(in->xsize,in->ysize);
  /* get memory for the image of gradient modulus */
  *grad = new_image_double(in->xsize,in->ysize);
  *gradx = new_image_double(in->xsize,in->ysize);
  *grady = new_image_double(in->xsize,in->ysize);

  /* get memory for "ordered" coordinate list */
  list = (struct coorlist *) calloc(xsize*ysize,sizeof(struct coorlist));
  *mem_p = (void *) list;
  range_l_s = (struct coorlist **) calloc(n_bins,sizeof(struct coorlist *));
  range_l_e = (struct coorlist **) calloc(n_bins,sizeof(struct coorlist *));
  if( list == NULL || range_l_s == NULL || range_l_e == NULL )
    error("list: not enough memory");
  for(i=0;i<n_bins;i++) range_l_s[i] = range_l_e[i] = NULL;

  /* 'undefined' on the down and right boundaries */
  for(i=0;i<ysize;i++) angles->data[i*xsize+xsize-1] = NOTDEF;
  for(j=0;j<xsize;j++) angles->data[(ysize-1)*xsize+j] = NOTDEF;

  /*** remaining part ***/
  for(i=0;i<ysize-1;i++)
    for(j=0;j<xsize-1;j++)
      {
        adr = i*xsize+j;

        /*
           Norm 2 computation using 2x2 pixel window:
             A B
             C D
           and
             com1 = D-A,  com2 = B-C.
           Then
             gx = B+D - (A+C)   horizontal difference
             gy = C+D - (A+B)   vertical difference
           com1 and com2 are just to avoid 2 additions.
         */
        com1 = in->data[adr+xsize+1] - in->data[adr];
        com2 = in->data[adr+1] - in->data[adr+xsize];
        gx = com1+com2;
        gy = com1-com2;
        norm2 = gx*gx+gy*gy;
        norm = sqrt( norm2 / 4.0 );

        (*grad)->data[adr] = norm;
	(*gradx)->data[adr] = gx/2;
	(*grady)->data[adr] = gy/2;
        if(norm <= threshold) /* norm too small, gradient not defined */
            angles->data[adr] = NOTDEF;
        else
	  {
            /* angle computation */
	    angles->data[adr] = atan2(gy,gx);
            /* store the point in the right bin according to its norm */
            ind = (unsigned int) (norm * (double) n_bins / max_grad);
            if(ind >= n_bins) ind = n_bins-1;
            if(range_l_e[ind] == NULL)
                range_l_s[ind] = range_l_e[ind] = list+list_count++;
            else
	      {
                range_l_e[ind]->next = list+list_count;
                range_l_e[ind] = list+list_count++;
              }
            range_l_e[ind]->x = (int) j;
            range_l_e[ind]->y = (int) i;
            range_l_e[ind]->next = NULL;
         }
      }

  /* Make the list of points "ordered" by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with higher gradient value.
   */
  for(i=n_bins-1; i>0 && range_l_s[i]==NULL; i--);
  start = range_l_s[i];
  end = range_l_e[i];
  if(start != NULL)
    for(i--;i>0; i--)
      if( range_l_s[i] != NULL )
        {
          end->next = range_l_s[i];
          end = range_l_e[i];
        }
  *list_p = start;

  /* free memory */
  free( (void *) range_l_s );
  free( (void *) range_l_e );
  return angles;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Is point (x,y) aligned to angle theta, up to precision 'prec'?
 */
int isaligned( int x, int y, image_double angles, double theta,
                      double prec )
{
  double a;

  /* check parameters */
  if( angles == NULL || angles->data == NULL )
    error("isaligned: invalid image 'angles'.");
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("isaligned: (x,y) out of the image.");
  if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");

  /* angle at pixel (x,y) */
  a = angles->data[ x + y * angles->xsize ];

  /* pixels whose level-line angle is not defined
     are considered as NON-aligned */
  if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
                                      'double_equal' here because there is
                                      no risk of problems related to the
                                      comparison doubles, we are only
                                      interested in the exact NOTDEF value */

  /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
  theta -= a;
  if( theta < 0.0 ) theta = -theta;
  if( theta > M_3_2_PI )
    {
      theta -= M_2__PI;
      if( theta < 0.0 ) theta = -theta;
    }

  return theta <= prec;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Group neighbour pixels sharing the same orientation up to precision 'prec'.
 */
void region_grow( int x, int y, image_double angles, struct point * reg,
                         int * reg_size, double * reg_angle, image_char used,
                         double prec )
{
  double sumdx,sumdy;
  int xx,yy,i;
  /* check parameters */
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("region_grow: (x,y) out of the image.");
  if( angles == NULL || angles->data == NULL )
    error("region_grow: invalid image 'angles'.");
  if( reg == NULL ) error("region_grow: invalid 'reg'.");
  if( reg_size == NULL ) error("region_grow: invalid pointer 'reg_size'.");
  if( reg_angle == NULL ) error("region_grow: invalid pointer 'reg_angle'.");
  if( used == NULL || used->data == NULL )
    error("region_grow: invalid image 'used'.");

  /* first point of the region */
  /* *reg_size = 1; */
  reg[*reg_size-1].x = x;
  reg[*reg_size-1].y = y;
  *reg_angle = angles->data[x+y*angles->xsize];  /* region's angle */
  sumdx = cos(*reg_angle);
  sumdy = sin(*reg_angle);
  used->data[x+y*used->xsize] = USED;
 
  /* try neighbors as new region points */
  for(i=*reg_size-1; i<*reg_size; i++)
    for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
      for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
        if( xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
            used->data[xx+yy*used->xsize] != USED &&
            isaligned(xx,yy,angles,*reg_angle,prec) )
          {
            /* add point */
            used->data[xx+yy*used->xsize] = USED;            
            reg[*reg_size].x = xx;
            reg[*reg_size].y = yy;
            ++(*reg_size);            
            
            /* update region's angle */
            sumdx += cos( angles->data[xx+yy*angles->xsize] );
            sumdy += sin( angles->data[xx+yy*angles->xsize] );
            *reg_angle = atan2(sumdy,sumdx);
          }
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** 
 */
void init_rect(struct rect *rec)
{
  rec->x1 = 0.0; rec->y1 = 0.0;
  rec->x2 = 0.0; rec->y2 = 0.0;
  rec->x  = 0.0; rec->y  = 0.0;
  rec->dx = 0.0; rec->dy = 0.0;
  rec->width = 0.0;
  rec->theta = 0.0;
  rec->p = 0.0; rec->prec = 0.0;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** 
 */
void EllipseDetection(image_double image,double rho,double prec,double p, 
                      double eps,int smooth,int *ell_count, 
                      int *circ_count,int *line_count,char *fstr)
{
  image_double angles,gradx,grady,grad,imgauss;
  image_char used;
  void *mem_p;
  int n_bins = 1024;
  double max_grad = 255.0;
  struct coorlist *list_p;
  struct point *reg, *regl;
  struct point3 *regc, *rege;
  struct rect rec;
  int reg_size = 0,regp_size[3];
  unsigned int xsize,ysize; /* image size */
  int i;
  int min_size[3];
  FILE *svg;
  double logNT[3]; /* number of tests for the 3 primitive types */ 
  double nfa[3]; /* NFA value using the discrete formulation for the three primitive types */
  double parame[5], paramc[5]; /* ellipse/circle parameters */
  double lin[5];
  double density_th = 0.7;
  double mlog10eps = - log10(eps); 
  double reg_angle;
  int pext[8];
  unsigned int xsz0,ysz0;
  char svgname[300];
  xsz0 = image->xsize;
  ysz0 = image->ysize;

  /* perform gaussian smoothing and subsampling */
  if (smooth)
    {
      imgauss = gaussian_sampler( image, 0.8, 0.6);
      free_image_double(image);
      /* compute gradient magnitude and orientation  */
      angles = ll_angle(imgauss,rho,&list_p,&mem_p,&gradx,&grady,&grad,n_bins,max_grad);
    }
  else
    angles = ll_angle(image,rho,&list_p,&mem_p,&gradx,&grady,&grad,n_bins,max_grad);

  xsize = angles->xsize;
  ysize = angles->ysize;

  /* display detection result */
  sprintf(svgname,"%s.svg",fstr);
  svg = init_svg(svgname,xsz0,ysz0);
  
  /* number of tests for elliptical arcs */
  logNT[2] = 4.0 *(log10((double)xsize)+log10((double)ysize)) + log10(9.0) + log10(3.0); /* N^8 */
  /* number of tests for circular arcs */
  logNT[1] = 3.0 *(log10((double)xsize)+log10((double)ysize)) + log10(9.0) + log10(3.0); /* N^6 */
  /* number of tests for line-segments */
  logNT[0] = 5.0 *(log10((double)xsize)+log10((double)ysize))/2.0 + log10(11) + log10(3.0); /* N^5 */

  /* thresholds from which an elliptical/circular/linear arc could be meaningful */
  min_size[2] =(int)((-logNT[2]+log10(eps))/log10(p));
  min_size[1] =(int)((-logNT[1]+log10(eps))/log10(p));
  min_size[0] =(int)((-logNT[0]+log10(eps))/log10(p));

  /* file to write coordinates of detected ellipses */
  FILE *fe = fopen("ellipses.txt","wt");

  /* allocate memory for region lists */
  reg = (struct point *) calloc(xsize * ysize, sizeof(struct point));
  regl = (struct point *) calloc(xsize * ysize, sizeof(struct point));
  regc = (struct point3 *) calloc(xsize * ysize, sizeof(struct point3));
  rege = (struct point3 *) calloc(xsize * ysize, sizeof(struct point3));
  used = new_image_char_ini(xsize,ysize,NOTUSED);
  
  /* init temporary buffers */
  gBufferDouble = (double*)malloc(sizeof(double));
  gBufferInt    = (int*)malloc(sizeof(int));

  init_rect(&rec);

  /* begin primitive detection */
  for(;list_p; list_p = list_p->next)
    { 
      reg_size = 0;
      if(used->data[list_p->y*used->xsize+list_p->x]==NOTUSED &&
         angles->data[list_p->y*angles->xsize+list_p->x] != NOTDEF)
	{
          /* init some variables */ 	
          for (i=0;i<5;i++) 
            {
              parame[i] = 0.0; paramc[i] = 0.0;
            }
          nfa[2] = nfa[1] = nfa[0] = mlog10eps; 
          reg_size = 1;regp_size[0] = regp_size[1] = regp_size[2] = 0;
	  region_grow(list_p->x, list_p->y, angles, reg, &reg_size, &reg_angle,
                          used, prec); 
                               

	  /*-------- FIT A LINEAR SEGMENT AND VERIFY IF VALID ------------------- */
	  valid_line(reg,&reg_size,reg_angle,prec,p,&rec,lin,grad,gradx,grady, 
	             used,angles,density_th,logNT[0],mlog10eps,&nfa[0]); 
          regp_size[0] = reg_size;

          for (i=0;i<regp_size[0];i++) {regl[i].x = reg[i].x; regl[i].y = reg[i].y; }

          if (reg_size>2)
            {
  
              /*-------- FIT CONVEX SHAPES (CIRCLE/ELLIPSE) AND VERIFY IF VALID -------- */
              valid_curve(reg,&reg_size,prec,p,angles,used,grad,gradx,grady,paramc,parame,
                          rec,logNT,mlog10eps,density_th,min_size,nfa,pext,regc,rege,regp_size);
              
              /* ------ DECIDE IF LINEAR SEGMENT OR CIRCLE OR ELLIPSE BY COMPARING THEIR NFAs -------*/
              if(nfa[2]>mlog10eps && nfa[2]>nfa[0] && nfa[2]>nfa[1] && regp_size[2]>min_size[2]) /* ellipse */
                {
                  (*ell_count)++;
                  /*if (smooth)
                    fprintf(fe,"%f %f %f %f %f \n",parame[0]*1.25,parame[1]*1.25,parame[2]*1.25,parame[3]*1.25,parame[4]);
                  else  
                    fprintf(fe,"%f %f %f %f %f \n",parame[0],parame[1],parame[2],parame[3],parame[4]);*/
                  write_svg_ellipse(fe,svg,parame,pext,smooth);
                  for (i=0;i<regp_size[0];i++)
                        used->data[regl[i].y*used->xsize+regl[i].x] = NOTUSED;
                  for (i=0;i<regp_size[1];i++)
                        used->data[regc[i].y*used->xsize+regc[i].x] = NOTUSED;
                  for (i=0;i<regp_size[2];i++)
                      if (rege[i].z == USEDELL)
                        { 
                        used->data[rege[i].y*used->xsize+rege[i].x] = USED;}
                      else{ 
                        used->data[rege[i].y*used->xsize+rege[i].x] = USEDELLNA;}
                }
              else if(nfa[1]>mlog10eps && nfa[1]>nfa[0] && nfa[1]>nfa[2] && regp_size[1]>min_size[1]) /* circle */
                     {
                       (*circ_count)++;
                       /*if (smooth)
                         fprintf(fe,"%f %f %f %f %f \n",paramc[0]*1.25,paramc[1]*1.25,paramc[2]*1.25,paramc[3]*1.25,paramc[4]);
                       else  
                         fprintf(fe,"%f %f %f %f %f \n",paramc[0],paramc[1],paramc[2],paramc[3],paramc[4]);*/
                       write_svg_circle(fe,svg,paramc,pext,smooth);
                       for (i=0;i<regp_size[0];i++)
                         used->data[regl[i].y*used->xsize+regl[i].x] = NOTUSED;
                       for (i=0;i<regp_size[2];i++)
                         used->data[rege[i].y*used->xsize+rege[i].x] = NOTUSED;
                       for (i=0;i<regp_size[1];i++)
                         if (regc[i].z == USEDCIRC)
                           used->data[regc[i].y*used->xsize+regc[i].x] = USED;
                         else 
                           used->data[regc[i].y*used->xsize+regc[i].x] = USEDCIRCNA;
                     }
                   else if(nfa[0]>mlog10eps && regp_size[0]>min_size[0] && nfa[0]>nfa[1] && nfa[0]>nfa[2]) /* line */
                          {
                            (*line_count)++;
			    write_svg_line(svg,lin,smooth);                    
		            for (i=0;i<regp_size[1];i++)
                              used->data[regc[i].y*used->xsize+regc[i].x] = NOTUSED;
                            for (i=0;i<regp_size[2];i++)
                              used->data[rege[i].y*used->xsize+rege[i].x] = NOTUSED;
                            for (i=0;i<regp_size[0];i++)
                              used->data[regl[i].y*used->xsize+regl[i].x] = USED;
                          }
                        else /* no feature */
                          { 
                            for (i=0;i<regp_size[1];i++)
                              used->data[regc[i].y*used->xsize+regc[i].x] = NOTUSED;
                            for (i=0;i<regp_size[2];i++)
                              used->data[rege[i].y*used->xsize+rege[i].x] = NOTUSED;
                            for (i=0;i<regp_size[0];i++)
                              used->data[regl[i].y*used->xsize+regl[i].x] = NOTUSED;
                          }
            }
        }/* IF USED */
    }/* FOR LIST */		    



  if (smooth) free_image_double(imgauss);
  else free_image_double(image);
  free_image_double(gradx); free_image_double(grady);
  free_image_double(grad); free_image_double(angles); 
  free_image_char(used);
  free(reg);free(regl); free(regc); free(rege); 
  free(gBufferDouble); free(gBufferInt); 
  free(mem_p);
  fclose(fe);
  fclose_svg(svg);
}
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/** Check if ellipse axes are positive and ensure big axis is first.
 */
int check_ellipse(double *param)
{
  double tmp;
  if (param[2]<=0||param[3]<=0) return 0;
  if (param[2]<param[3])
    {
      tmp = param[3]; param[3] = param[2]; param[2] = tmp; param[4] += M_1_2_PI;
    }
  if (param[2]/param[3]>100) return 0;
  if (param[4]>M_1_2_PI) 
    param[4] = param[4] - M_PI;
  if (param[4]>M_PI) 
    param[4] = param[4] - M_PI;
  if (param[4]<-M_1_2_PI)
    param[4] = param[4] + M_PI;
  return 1;
}

/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  if (argc<2) error("use : ./elsd image_name.pgm");
  double quant = 2.0;       /* Bound to the quantization error on the
                                gradient norm.                                */
  double ang_th = 22.5;     /* Gradient angle tolerance in degrees.           */
  double p = ang_th/180.0;
  double prec = M_PI*ang_th/180.0; /* radian precision */
  double rho = quant/sin(prec);
  double eps = 1; //atof(argv[2]);
  int smooth = 1; //atoi(argv[3]);
  int ell_count = 0, line_count = 0, circ_count = 0;
  image_double image;

  image = read_pgm_image_double(argv[1]);

  EllipseDetection(image, rho, prec, p, eps, smooth, &ell_count, &circ_count, 
                   &line_count,argv[1]);
  printf("%s\n", argv[1]);
  printf("%d elliptical arcs, %d circular arcs, %d line segments\n", 
         ell_count, circ_count, line_count);
  return 0;
}
/*----------------------------------------------------------------------------*/


