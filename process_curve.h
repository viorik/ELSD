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

#ifndef PROCESS_CURVE_H
#define PROCESS_CURVE_H

void valid_curve(struct point *reg,int *reg_size,double prec,double p,
                 image_double angles,image_char used,image_double grad,
                 image_double gradx,image_double grady,double *paramc,
                 double *parame,struct rect rec,double *logNT,double mlog10eps,
                 double density_th,int *min_size,double *nfa,int *pext,
                 struct point3 *regc,struct point3 *rege,int *regp_size);
double d_rosin (double *param,double x,double y);
double angle(double x,double y,double *foci);
void param2ellipse(double *param, double *ellipse);
void rosin_point (double *param,double x,double y,double *xi,double *yi);
void fitellipse(int reg_size,double *vgg,double *param);
void fitcircle(int reg_size,double *vgg,double *param);
void fit_equations(struct point *reg,int reg_size,image_double gradx, 
                   image_double grady,double *vgg);
#endif
