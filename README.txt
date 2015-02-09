ABOUT THIS SOURCE CODE
The files in this folder contain the source code of ELSD, published in 
'A Parameterless Line Segment and Elliptical Arc Detector with Enhanced Ellipse
Fitting', V. Patraucean, P. Gurdjos, R. Grompone von Gioi, ECCV2012.

Corresponding author: viorica patraucean vpatrauc@enseeiht.fr.  
 
The code generating and validating line segment hypotheses is taken up from 
LSD source code, available at http://www.ipol.im/pub/art/2012/gjmr-lsd/.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any
later version. 


REQUIREMENTS
The ELSD source code needs the CLAPACK/CBLAS library for some linear algebra 
computations. Version 3.2.1 was used.


SOURCE CODE FILES
'elsd.c'          contains the main() function, some general functions (e.g. 
                  read/write pgm images), and functions shared by the line 
                  segment detection part and the ellipse detection part (e.g. 
                  NFA computation);
'process_line.c'  contains functions used to produce and validate line segment 
                  hypotheses; the code is taken up mainly from LSD source code.
'process_curve.c' contains functions used to produce ellipse/circle hypotheses;
'valid_curve.c'   contains functions to validate circle/ellipse hypotheses;
'write_svg.c'     functions to write the result in svg format.


COMPILATION
'makefile' 	  example of makefile to compile the source code. If the paths 
                  to the libraries are ok, a simple 'make' would compile the 
                  code and produce the executable called 'elsd'.


EXECUTION
./elsd imagename  runs ELSD on the image specified by 'imagename'. This ELSD 
                  version works only with PGM images. This folder contains the
                  image 'stars.pgm' for testing purposes.   


OUTPUT
'imagename.svg'   contains the execution result in SVG format. 'stars.pgm.svg' 
                  contains the result for the sample image 'stars.pgm'.   	  
'ellipses.txt'    contains the parameters of the detected circular/elliptical 
                  arcs in the form 'x_c y_c a b theta ang_start ang_end'.
In the console, the name of the output svg file is displayed, and the number
of features of each type. For the 'stars.pgm' image, the output is
stars.pgm.svg
17 elliptical arcs, 310 circular arcs, 165 line segments.

 
