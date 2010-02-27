/*   qrs/main.c                                                            *
 *   this does a cubic spline approximation.                               *
 ***************************************************************************
 *   Copyright (C) 2009, 2010 Hermann Sonntag                              *
 *   hermann.sonntag@tu-ilmenau.de                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cspl/cspl.h>

//#include "config.h"
#include "cmdline.h"

/* constants */

/* options */

#define ARG_NONE 0
#define ARG_REQUIRED 1
#define ARG_OPTIONAL 2

int read_signals(void * dat, size_t size, size_t channels, size_t count, char *filename) {
    FILE *fp;

    /* open the file */
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        perror ("The following error occurred in read_signals");
        return 1;
    }

    /* copy the file into memory */
    if (fread(dat, size, channels*count, fp) != channels*count){
        perror ("The following error occurred in read_signals");
        return 2;
    }

    /* close the file */
    return fclose (fp);
}

/*!
 * This is the main method of the QRS Splineapproximation.
 *
 * @author Hermann Sonntag
 *
 */
int main(int argc, char *argv[]) {

    struct gengetopt_args_info args_info;
    char * filename = NULL;
    unsigned int i;

    /* let's call our cmdline parser */
    if (cmdline_parser (argc, argv, &args_info) != 0)
        exit(1) ;

    for ( i = 0 ; i < args_info.inputs_num ; ++i ) {
        filename = (char *) malloc (strlen (args_info.inputs[i]));
        strcpy(filename, args_info.inputs[i]);
    }

    size_t splinelength = 0;
    if ( args_info.int_opt_given )
        splinelength = args_info.int_opt_arg;

    cmdline_parser_free (&args_info);
    double t[splinelength], x[splinelength], c_xy[splinelength];
    double qrs_tmpl[splinelength];
    unsigned int m_corr[splinelength];
    int interval;
    int _dat[2][splinelength];
    read_signals(_dat, 4, 2, splinelength, filename);
    for (i = 0; i < splinelength; i++) {
        x[i] = (double) _dat[1][i];
        t[i] = 0.001 * i;
        if (i < 100)
            qrs_tmpl[i] =_dat[1][i];
        else
            qrs_tmpl[i] = qrs_tmpl[99];
    }

    cspl_norm (x, splinelength);
    cspl_norm (qrs_tmpl, splinelength);
    /*for (i = 0; i < splinelength; i++) {
        printf("%d %f\n", i, x[i]);
    }*/

    int a = 0;
    if (args_info.enum_opt_given) {
        if (strcmp (args_info.enum_opt_arg, cmdline_parser_enum_opt_values[0]) == 0)
        {
            cspl_radix2_xcorr (c_xy, qrs_tmpl, x, splinelength);
            cspl_norm (c_xy, splinelength);
            for (i = 0; i < splinelength; i++) {
                //printf("%d %f\n", i, c_xy[i]);
            }
            a = cspl_eval_periodic_max (m_corr, c_xy, splinelength, 0.8);
        }
    }
    else
    {
        cspl_root_mean_square (c_xy, qrs_tmpl, x, 100, splinelength);
        cspl_norm (c_xy, splinelength);
        for (i = 0; i < splinelength; i++) {
            //printf("%d %f\n", i, c_xy[i]);
        }
        a = cspl_eval_periodic_min2 (m_corr, c_xy, splinelength, 200, -0.2);
    }

    //printf ("a=%d\n", a);
    /*printf ("#m=0,S=2\n");
    for (i = 0; i < a; i++) {
        printf ("%d %d\n",i, m_corr[i]);
    }*/
    for (i = 0; i < splinelength; i++)
        qrs_tmpl[i] = 0.0;

    interval = cspl_average (qrs_tmpl, x, m_corr, splinelength, a);
    cspl_norm (qrs_tmpl, interval);
    /*printf ("#m=1,S=0\n");
    for (i = 0; i < interval; i++) {
        printf("%.5f %.5f\n", t[i], qrs_tmpl[i]);
    }
    */
    cspl_qrs_init();

    gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, interval);
    gsl_spline_init(spline, t, qrs_tmpl, interval);

    double sigma[interval];
    const size_t p = 4;
    const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc (T, interval, p);
    gsl_matrix * covar = gsl_matrix_alloc (p, p);
    double x_init[4] = { 1.0, 0.0, 0.0, 0.0 };
    gsl_vector_view vec = gsl_vector_view_array (x_init, p);
    double params[p][a];
    for (i = 0; i < a; i++) {
        struct cspl_qrs_data data = {interval, t, x + m_corr[i], spline, sigma, p, s, covar, vec};
        cspl_qrs_fit (&data);
        //size_t j;
        params[0][i] = gsl_vector_get(s->x, 0);
        params[1][i] = gsl_vector_get(s->x, 1);
        params[2][i] = gsl_vector_get(s->x, 2);
        params[3][i] = gsl_vector_get(s->x, 3);

        /*printf ("#m=0,S=0\n");
          for (j = 0; j < interval; j++) {
          double y_j = gsl_vector_get(s->x, 0)*gsl_spline_eval(spline, t[j] - gsl_vector_get(s->x, 1), acc) + gsl_vector_get(s->x, 2)*t[j] + gsl_vector_get(s->x, 3);
          printf ("%.5f %.5f\n", t[j + m_corr[i]], y_j);
          }
          printf ("#m=3,S=0\n");
          for (j = m_corr[i]; j < m_corr[i] + interval; j++) {
          printf ("%.5f %.5f\n", t[j], x[j]);
          }
          */
    }
    cspl_real_fft (params[0], a);
    cspl_real_fft (params[1], a);
    cspl_real_fft (params[2], a);
    cspl_real_fft (params[3], a);
    printf ("#m=0,S=0\n");
    for (i = 0; i < a; i++) {
        printf("%.5f %.5f\n", t[m_corr[i]], params[0][i]);
    }
    printf ("#m=1,S=0\n");
    for (i = 0; i < a; i++) {
        printf("%.5f %.5f\n", t[m_corr[i]], params[1][i]);
    }
    printf ("#m=2,S=0\n");
    for (i = 0; i < a; i++) {
        printf("%.5f %.5f\n", t[m_corr[i]], params[2][i]);
    }
    printf ("#m=3,S=0\n");
    for (i = 0; i < a; i++) {
        printf("%.5f %.5f\n", t[m_corr[i]], params[3][i]);
    }
    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    cspl_qrs_free ();
    return EXIT_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f"
            "|f(x)| = %g\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            gsl_blas_dnrm2 (s->f));
}

