/*     splineapprox does a cubic spline approximation
 *     Copyright (C) 2009 Hermann Sonntag
  
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.

 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.

 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>

//#include "config.h"
#include "cmdline.h"

#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "qrsfit.c"
     
#define N 128
     
/* constants */

#define SPLINELENGTH 128
#define CHANNELS 8
#define DATASIZE 2
#define SAMPLING_TIME 0.1
#define T_QRS 3.1415

/* options */

#define	ARG_NONE	0
#define	ARG_REQUIRED	1
#define	ARG_OPTIONAL	2

void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int read_signals(double ***dat, unsigned short size, unsigned short channels, unsigned int count, char *filename) {
  FILE *fp;

  /* open the file */
  fp = fopen(filename, "rb");
  if (fp == NULL) {
    perror ("The following error occurred in read_signals");
    return(1);
  }
  
  /* copy the file into memory */
  int i;
  for (i = 0; i < count; i++) {
    if(fread(dat[i], size, channels, fp) != channels){
      perror ("The following error occurred in read_signals");
      return(2);
    }
  }

  /* close the file */
  return fclose(fp);
}

/*!
 * This is the main method of the QRS Splineapproximation.
 *
 * @author Hermann Sonntag
 *
 */
int main(int argc, char *argv[]) {

  struct gengetopt_args_info args_info;
  char filename[64];
  unsigned int i;
     
  /*        printf( "This one is from a C program \n" ); */
  /*        printf( "Try to launch me with some options\n" ); */
  /*        printf( "(type splineapprox --help for the complete list)\n" ); */
  /*        printf( "For example: ./splineapprox *.* --funct-opt\n" ); */
       
  /*        /\* let's call our cmdline parser *\/ */
  /*        if (cmdline_parser (argc, argv, &args_info) != 0) */
  /*          exit(1) ; */
     
  /*        printf( "Here are the options you passed...\n" ); */
     
  /*        for ( i = 0 ; i < args_info.inputs_num ; ++i ) { */
  /*          printf( "file: %s", args_info.inputs[i] ); */
  /*          strcpy(filename, args_info.inputs[i]); */
  /*        } */

  /*        if ( args_info.funct_opt_given ) */
  /*          printf("You chose --funct-opt or -F." ); */
     
  /*        if ( args_info.str_opt_given ) */
  /*          printf( "You inserted %s%s%s", args_info.str_opt_arg, " for ", "--str-opt option." ); */
     
  /*        if ( args_info.int_opt_given ) */
  /*          printf( "This is the integer you input: %d%s ", args_info.int_opt_arg, "." ); */
     
  /*        if (args_info.flag_opt_given) */
  /*          printf( "The flag option was given! " ); */
     
  /*        printf( "The flag is %s%s", ( args_info.flag_opt_flag ? "on" : "off" ), ". " ); */
     
  /*        if (args_info.enum_opt_given) { */
  /*          printf( "enum-opt value: %s", args_info.enum_opt_arg ); */
  /*          printf( "enum-opt (original specified) value: %s", args_info.enum_opt_orig ); */
  /*        } */
     
  /*        if (args_info.secret_given) */
  /*          printf( "Secret option was specified: %d", args_info.secret_arg ); */
     
  /*        printf( "%s! ", args_info.def_opt_arg ); */
     
  /*        printf( "Have a nice day! :-)\n" ); */
     
  /*        cmdline_parser_free (&args_info); /\* release allocated memory *\/ */


  double t[SPLINELENGTH], xs[SPLINELENGTH], y[SPLINELENGTH], c_xy[SPLINELENGTH], n_QRS[SPLINELENGTH] = {0}; 
  unsigned int m_corr[SPLINELENGTH] = {0}; 

  for (i = 0; i < SPLINELENGTH; i++) {
    t[i] = i * SAMPLING_TIME;
    xs[i] = sin(t[i]);
    if (t[i] < T_QRS)
      y[i] = xs[i];
    else
      y[i] = 0;
  }
 

  return EXIT_SUCCESS;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
	  "|f(x)| = %g\n",
	  iter,
	  gsl_vector_get (s->x, 0),
	  gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
	  gsl_blas_dnrm2 (s->f));
}
