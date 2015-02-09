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

#ifndef VALID_CURVE_H
#define VALID_CURVE_H

double valid_ellipse(struct point *reg, int reg_size, image_char used, 
                     double prec, double p, image_double angles,
                     image_double grad,image_double gradx,image_double grady,
                     double *param, double logNTE,int dir,
                     int *pext,struct point3 *rege, int *rege_size,int min_size,double mlog10eps);
double valid_circle(struct point *reg, int reg_size, image_char used, 
                    double prec, double p, image_double angles, 
                    image_double grad,image_double gradx,image_double grady,
                    double *param, double logNTC, int dir, 
                    int *pext,struct point3 *regc, int *regc_size,int min_size,double mlog10eps);
int isInAng(double ang, double ang1, double ang2);

#endif
