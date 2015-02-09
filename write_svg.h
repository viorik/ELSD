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

#ifndef WRITE_SVG_H
#define WRITE_SVG_H

FILE* init_svg(char * filename, unsigned int xsize, unsigned int ysize );
void fclose_svg(FILE *svg);
void write_svg_ellipse(FILE *fe, FILE *svg,double *param,int *pext,int smooth);
void write_svg_circle(FILE *fe, FILE *svg,double *param,int *pext,int smooth);
void write_svg_line(FILE *svg, double *lin,int smooth);

#endif
