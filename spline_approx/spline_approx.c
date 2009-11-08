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

#include "splineapprox.h"
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
  int z;
  z = gsl_fft_real_radix2_transform(xs, 1, SPLINELENGTH);
  z = gsl_fft_real_radix2_transform(y, 1, SPLINELENGTH);

  size_t n = SPLINELENGTH;
  c_xy[0] = xs[0]*y[0];
  c_xy[n/2] = xs[n/2]*y[n/2];
  for (i = 1; i < n/2; ++i) {
    c_xy[i] = xs[i]*y[i] + xs[n-i]*y[n-i];
    c_xy[n-i] = xs[n-i]*y[i] - xs[i]*y[n-i];
  }
  z = gsl_fft_halfcomplex_radix2_inverse(c_xy, 1, n);

  for (i = 0; i < SPLINELENGTH; i++) {
    t[i] = i * SAMPLING_TIME;
    xs[i] = sin(t[i]);
    if (t[i] < T_QRS)
      y[i] = xs[i];
    else
      y[i] = 0;
  }

  size_t m;
  double max;

  for (i = 0; i < n; i++) {
    max = 0;
    for (m = m_corr[i]; (m < (m_corr[i] + (int)(1.5*T_QRS/SAMPLING_TIME))) && (m < n); m++) {
      if (c_xy[m] > max) {
	max = c_xy[m];
	m_corr[i] = m;
      }
    }
    if((m == n) || (i == (n - 2)))
      break;
    m_corr[i + 1] = m_corr[i] +  (int)(0.5*T_QRS/SAMPLING_TIME);
  }

  max = 0;
  size_t count = i;
  for (i = 0; i < count; i++) {
    for (m = m_corr[i]; m < m_corr[i] + T_QRS/SAMPLING_TIME; m++) {
      n_QRS[m - m_corr[i]] += xs[m];
      if (n_QRS[m - m_corr[i]]>max)
	max = n_QRS[m - m_corr[i]];
    }
  }

  for (i = 0; i < n; i++) {
    n_QRS[i] /= max;
  }

  /*------------------------*/

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  //const size_t n = N;
  const size_t p = 3;
     
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double sigma[N];
  struct data d = { n, 0, t, y, n_QRS, sigma};
  gsl_multifit_function_fdf f;
  double x_init[3] = { 2.0, 0.0, 2.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;
     
  gsl_rng_env_setup();
     
  type = gsl_rng_default;
  r = gsl_rng_alloc (type);
     
  f.f = &s_QRS_f;
  f.df = &s_QRS_df;
  f.fdf = &s_QRS_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;
     
  /* This is the data to be fitted */
     
  for (i = 0; i < n; i++)
    {
      /*  double t = i; */
      /*       y[i] = 1.0 + 5 * exp (-0.1 * t)  */
      /* 	+ gsl_ran_gaussian (r, 0.1); */
      sigma[i] = 0.1;
      /*  printf ("data: %u %g %g\n", i, y[i], sigma[i]); */
    };
     
  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
     
  print_state (iter, s);
     
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
     
      printf ("status = %s\n", gsl_strerror (status));
     
      print_state (iter, s);
     
      if (status)
	break;
     
      status = gsl_multifit_test_delta (s->dx, s->x,
					1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);
     
  gsl_multifit_covar (s->J, 0.0, covar);
     
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
     
  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
     
    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof); 
     
    printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0)); 
    printf ("t_beat = %.5f +/- %.5f\n", FIT(1), c*ERR(1)); 
    printf ("S_0      = %.5f +/- %.5f\n", FIT(2), c*ERR(2)); 
  }
     
  printf ("status = %s\n", gsl_strerror (status));
     
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);


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
