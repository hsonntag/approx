/*   qrs/main.c                                                            *
 *   this does a cubic spline approximation and analysis.                  *
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
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <cspl/cspl.h>

#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

//#include "config.h"
#include "cmdline.h"

/* constants */

/* options */

#define ARG_NONE 0
#define ARG_REQUIRED 1
#define ARG_OPTIONAL 2

int read_signals(int * dat, size_t size, size_t channels, size_t count, char *filename) {
    FILE *fp;

    /* open the file */
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        perror ("The following error occurred in read_signals");
        return 1;
    }

    /* copy the file into memory */
    int i;
    for (i = 0; (i < count) && (!feof(fp)); i++)
    {
        if (fread(&dat[i], size, 1, fp) != 1)
        {
            perror ("The following error occurred in read_signals");
            return 2;
        }
        fseek(fp, size*(206-1), SEEK_CUR);
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
    int timedomain = 0;
    int template = 0;
    uint8_t time = 0;
    uint8_t frq = 0;

    struct gengetopt_args_info args_info;
    char * filename = NULL;
    unsigned int i, j;

    /* let's call our cmdline parser */
    if (cmdline_parser (argc, argv, &args_info) != 0)
        exit(1) ;

    for ( i = 0 ; i < args_info.inputs_num ; ++i ) {
        filename = (char *) malloc (strlen (args_info.inputs[i]));
        strcpy (filename, args_info.inputs[i]);
    }

    size_t splinelength = 0;
    if ( args_info.int_opt_given )
    {
        splinelength = args_info.int_opt_arg;
    }

    if ( args_info.frq_opt_given )
    {
        frq = args_info.frq_opt_arg;
    }

    if ( args_info.time_opt_given )
    {
        time = args_info.time_opt_arg;
    }

    cmdline_parser_free (&args_info);
    double t[splinelength], x[splinelength], c_xy[splinelength];
    double * qrs_tmpl;
    qrs_tmpl = (double *) malloc (100*sizeof (double)); //TODO dynamic
    unsigned int * m_corr, * m_corr1, * m_corr2;
    int interval;
    int _dat[splinelength];
    read_signals (_dat, 4, 206, splinelength, filename);
    for (i = 0; i < splinelength; i++) {
        x[i] = (double) _dat[i];
        t[i] = ((double)1)/1025 * i;
        if (i < 100)
            qrs_tmpl[i] =_dat[i];
    }

    cspl_norm (x, splinelength);
    cspl_norm (qrs_tmpl, 100);

    //for (i = 0; i < splinelength; i++) {
    //    printf("%d %f\n", i, x[i]);
    //}


    int a = 0;
    int b = 0;
    int k = 0;
    //if (args_info.enum_opt_given)
    {
        //  if (strcmp (args_info.enum_opt_arg, cmdline_parser_enum_opt_values[0]) == 0)
        {
            cspl_xcorr (c_xy, qrs_tmpl, x, 100, splinelength);
            cspl_norm (c_xy, splinelength);
            /*
               for (i = 0; i < splinelength; i++) {
               printf("%d %f\n", i, c_xy[i]);
               }
               */
            m_corr1 = (unsigned int *) malloc (sizeof (unsigned int) * splinelength / 200);
            a = cspl_eval_periodic_max2 (m_corr1, c_xy, splinelength, 200, 0.15);
        }
    }
    //else
    {
        cspl_root_mean_square (c_xy, qrs_tmpl, x, 100, splinelength);
        cspl_norm (c_xy, splinelength);
        /* 
           for (i = 0; i < splinelength; i++) {
           printf("%d %f\n", i, c_xy[i]);
           }
           */
        m_corr2 = (unsigned int *) malloc (sizeof (unsigned int) * splinelength / 200);
        b = cspl_eval_periodic_min2 (m_corr2, c_xy, splinelength, 200, -0.15);
    }
    m_corr = (unsigned int *) malloc (sizeof (unsigned int) * GSL_MIN (a, b));
    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++)
        {
            if (m_corr1[i] == m_corr2[j])
            {
                m_corr[k] = m_corr1[i];
                k++;
            }
        }
    }

    printf ("#m=0,S=2\n");
    for (i = 0; i < k; i++) 
    {
        printf ("%d, %d\n", m_corr1[i], m_corr2[i]);
        printf ("%d %d\n",i, m_corr[i]);
    }
    if (k != 0) 
    {
        qrs_tmpl = realloc (qrs_tmpl, splinelength/k*sizeof (double));
        for (i = 0; i < splinelength/k; i++)
        {
            qrs_tmpl[i] = 0.0;
        }
    }
    interval = cspl_average (qrs_tmpl, x, m_corr, splinelength, k);
    cspl_norm (qrs_tmpl, interval);
    if (template)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < interval; i++) {
            printf("%.5f %.5f\n", t[i], qrs_tmpl[i]);
        }
    }
    free (m_corr);
    if (a > b)
    {
        k = a;
        m_corr = m_corr1;
    }
    else
    {
        k = b;
        m_corr = m_corr2;
    }
    cspl_qrs_init();

    gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, interval);
    gsl_spline_init(spline, t, qrs_tmpl, interval);
    if (template) {
        printf ("#m=1,S=0\n");
        for (j = 0; j < interval*10; j++) 
            printf ("%.5f %.5f\n", t[j]/10, gsl_spline_eval(spline, t[j]/10, acc));
    }

    double sigma[interval];
    const size_t p = 4;
    const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc (T, interval, p);
    gsl_matrix * covar = gsl_matrix_alloc (p, p);
    double x_init[4] = { 1.0, 0.0, 0.0, 0.0 };
    gsl_vector_view vec = gsl_vector_view_array (x_init, p);
    double params[p][k];

    nfft_plan pl[p];

    /** init an one dimensional plan */
    nfft_init_1d(&pl[0],k,k);
    nfft_init_1d(&pl[1],k,k);
    nfft_init_1d(&pl[2],k,k);
    nfft_init_1d(&pl[3],k,k);
    for (i = 0; i < k; i++) {
        pl[0].x[i] = ((double) m_corr[i])/(splinelength);//-0.5;
        pl[1].x[i] = ((double) m_corr[i])/(splinelength);//-0.5;
        pl[2].x[i] = ((double) m_corr[i])/(splinelength);//-0.5;
        pl[3].x[i] = ((double) m_corr[i])/(splinelength);//-0.5;
    }
    /** precompute psi, the entries of the matrix B */
    //if(pl.nfft_flags & PRE_ONE_PSI)
    //    nfft_precompute_one_psi(&pl);
    for (i = 0; i < k; i++) {
        struct cspl_qrs_data data = {interval, t, x + m_corr[i], spline, sigma, p, s, covar, vec};
        cspl_qrs_fit (&data);
        size_t j;
        params[0][i] = gsl_vector_get(s->x, 0);
        params[1][i] = gsl_vector_get(s->x, 1);
        params[2][i] = gsl_vector_get(s->x, 2);
        params[3][i] = gsl_vector_get(s->x, 3);
        pl[0].f[i] = gsl_vector_get(s->x, 0);
        pl[1].f[i] = gsl_vector_get(s->x, 1);
        pl[2].f[i] = gsl_vector_get(s->x, 2);
        pl[3].f[i] = gsl_vector_get(s->x, 3);

        if (timedomain)
            printf ("#m=1,S=0\n");
        for (j = 0; j < interval; j++) 
        {
            double y_j = gsl_vector_get(s->x, 0)*gsl_spline_eval(spline, t[j] - gsl_vector_get(s->x, 1), acc) + gsl_vector_get(s->x, 2)*t[j] + gsl_vector_get(s->x, 3);
            if (timedomain)
                printf ("%.5f %.5f\n", t[j], y_j);
        }
        if (timedomain)
            printf ("#m=1,S=0\n");
        for (j = m_corr[i]; j < m_corr[i] + interval; j++) 
        {
            if (timedomain)
                printf ("%.5f %.5f\n", t[j], x[j]);
        }
    }

    /*
       cspl_real_fft (params[0], k);
       cspl_real_fft (params[1], k);
       cspl_real_fft (params[2], k);
       cspl_real_fft (params[3], k);
       cspl_halfcomplex_abs (params[0], k);
       cspl_halfcomplex_abs (params[1], k);
       cspl_halfcomplex_abs (params[2], k);
       cspl_halfcomplex_abs (params[3], k);
       */

    ndft_adjoint(&pl[0]);
    ndft_adjoint(&pl[1]);
    ndft_adjoint(&pl[2]);
    ndft_adjoint(&pl[3]);

    if (frq & 0x8)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1] - k/(2*t[splinelength - 1]), cabs(pl[0].f_hat[i]));
        }
    }
    if (frq & 0x4)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1] - k/(2*t[splinelength - 1]), cabs(pl[1].f_hat[i]));
        }
    }
    if (frq & 0x2)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1] - k/(2*t[splinelength - 1]), cabs(pl[2].f_hat[i]));
        }
    }
    if (frq & 0x1)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1] - k/(2*t[splinelength - 1]), cabs(pl[3].f_hat[i]));
        }
    }

    if (time & 0x8)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", t[m_corr[i]], params[0][i]);
        }
    }
    if (time & 0x4)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", t[m_corr[i]], params[1][i]);
        }
    }
    if (time & 0x2)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", t[m_corr[i]], params[2][i]);
        }
    }
    if (time & 0x1)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f\n", t[m_corr[i]], params[3][i]);
        }
    }
    /** finalise the one dimensional plan */
    nfft_finalize(&pl[0]);
    nfft_finalize(&pl[1]);
    nfft_finalize(&pl[2]);
    nfft_finalize(&pl[3]);

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);

    cspl_qrs_free ();

    free (m_corr1);
    free (m_corr2);
    free (qrs_tmpl);
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

