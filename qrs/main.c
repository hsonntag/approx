/*     Splineapprox does a cubic spline approximation.
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
#include <string.h>
#include <cspl/cspl.h>

//#include "config.h"
#include "cmdline.h"

/* constants */

/* options */

#define	ARG_NONE	0
#define	ARG_REQUIRED	1
#define	ARG_OPTIONAL	2

//void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int read_signals(void * dat, size_t size, size_t channels, size_t count, char *filename) {
	FILE *fp;

	/* open the file */
	fp = fopen(filename, "rb");
	if (fp == NULL) {
		perror ("The following error occurred in read_signals");
		return 1;
	}

	/* copy the file into memory */
	if(fread(dat, size, channels*count, fp) != channels*count){
		perror ("The following error occurred in read_signals");
		return 2;
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
	char * filename;
	unsigned int i;

	/*        printf( "This one is from a C program \n" ); */
	/*        printf( "Try to launch me with some options\n" ); */
	/*        printf( "(type splineapprox --help for the complete list)\n" ); */
	/*        printf( "For example: ./splineapprox *.* --funct-opt\n" ); */

	/*        /\* let's call our cmdline parser *\/ */
	if (cmdline_parser (argc, argv, &args_info) != 0) 
		exit(1) ; 

	/*        printf( "Here are the options you passed...\n" ); */

	for ( i = 0 ; i < args_info.inputs_num ; ++i ) { 
		filename = (char *) malloc(strlen(args_info.inputs[i]));
		strcpy(filename, args_info.inputs[i]); 
	} 

	/*        if ( args_info.funct_opt_given ) */
	/*          printf("You chose --funct-opt or -F." ); */

	/*        if ( args_info.str_opt_given ) */
	/*          printf( "You inserted %s%s%s", args_info.str_opt_arg, " for ", "--str-opt option." ); */
	size_t splinelength;
	if ( args_info.int_opt_given ) 
		splinelength = args_info.int_opt_arg;
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

	cmdline_parser_free (&args_info);
	double t[splinelength], x[splinelength], c_xy[splinelength]; 
	unsigned int m_corr[splinelength]; 
	int interval;
	int _dat[1][splinelength];
	read_signals(_dat, 4, 1, splinelength, filename);
	double qrs_tmpl[splinelength];
	for (i = 0; i < splinelength; i++) {
		x[i] = (double) _dat[0][i];
		t[i] = 1.0 * i;
		if (i < 100)
			qrs_tmpl[i] =_dat[0][i];
		else
			qrs_tmpl[i] = qrs_tmpl[99];
	}

	cspl_norm (x, splinelength);
	cspl_norm (qrs_tmpl, splinelength);
//	cspl_radix2_xcorr (c_xy, x, qrs_tmpl, splinelength);        
//	cspl_norm (c_xy, splinelength);
	cspl_root_mean_square (c_xy, qrs_tmpl, x, 100, splinelength);
	cspl_norm (c_xy, splinelength - 100);
	for (i = 0; i < splinelength; i++) {
		printf("%d %f\n", i, c_xy[i]);
	}
	int a =	cspl_eval_periodic_min2 (m_corr, c_xy, splinelength, 200, 0.09); 
		printf("a=%d\n", a);
	for (i = 0; i < a; i++) {
		printf("%d %d\n",i, m_corr[i]);
	}
	for (i = 0; i < splinelength; i++)
		qrs_tmpl[i] = 0.0;

interval = cspl_norm_average (qrs_tmpl, x, m_corr, splinelength, a);
		for (i = 0; i < splinelength; i++) {
			printf("%d %f\n", i, x[i]);
		}
		cspl_qrs_init();
		gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, interval);
		gsl_spline_init(spline, t, qrs_tmpl, interval);
		double sigma[interval];
//cspl_qrs_fit(t, x, spline, sigma, interval);
for (i = 0; i < a; i++)
cspl_qrs_fit_at(t, x, spline, sigma, interval, m_corr[i]);
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
