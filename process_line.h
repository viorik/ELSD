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

#ifndef PROCESS_LINE_H
#define PROCESS_LINE_H

void valid_line(struct point* reg,int *reg_size,double reg_angle,double prec,
                double p,struct rect *rec,double *l,image_double grad,
                image_double gradx,image_double grady, 
                image_char used,image_double angles,double density_th, 
                double logNT,double mlog10eps,double *nfaline);
int refine(struct point *reg,int reg_size0,int *reg_size,image_double modgrad,
           image_double gradx,image_double grady,double prec,double p,
           struct rect *rec,image_char used,image_double angles,
           double density_th,double logNT,double eps);
void region2rect(struct point *reg,int reg_size0,int reg_size,image_double modgrad,
                 double reg_angle, double prec, double p,struct rect * rec);
#endif

