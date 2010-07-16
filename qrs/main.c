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
#include <complex.h>
#include <cspl/cspl.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>

#include "nfft3util.h"
#include "nfft3.h"

//#include "config.h"
#include "cmdline.h"

/* constants */

/* options */

#define ARG_NONE 0
#define ARG_REQUIRED 1
#define ARG_OPTIONAL 2

int read_signals(void * dat, size_t size, size_t channel, size_t count, char *filename) {
FILE *fp;

/* open the file */
fp = fopen(filename, "rb");
if (fp == NULL) {
    perror ("The following error occurred in read_signals");
        return 1;
    }

    fseek(fp, size*(channel - 1), SEEK_CUR);
    /* copy the file into memory */
    int i;
    int ptr = (int)dat;
    for (i = 0; (i < count) && (!feof(fp)); i++)
    {
        if (fread((void *)ptr, size, 1, fp) != 1)
        {
            perror ("The following error occurred in read_signals");
            return 2;
        }
        ptr += size;
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
    int8_t timedomain = 0;
    int8_t errors = 0;
    int8_t template = 0;
    uint8_t time = 0;
    uint8_t frq = 0;

    unsigned int a, b, channel_start = 1, channel_count = 1,
                 k = 0, k_min = UINT_MAX, qrs_start = 0, samplingrate = 1025;
    size_t qrs_length = 100;
    struct gengetopt_args_info args_info;
    char * filename = NULL;
    unsigned int i, j, n;

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

    if ( args_info.qrs_start_opt_given )
    {
        qrs_start = args_info.qrs_start_opt_arg;
    }

    if ( args_info.qrs_length_opt_given )
    {
        qrs_length = args_info.qrs_length_opt_arg;
    }

    if ( args_info.channel_start_opt_given )
    {
        channel_start = args_info.channel_start_opt_arg;
    }

    if ( args_info.channel_count_opt_given )
    {
        channel_count = args_info.channel_count_opt_arg;
    }

    if ( args_info.samplingrate_opt_given )
    {
        samplingrate = args_info.samplingrate_opt_arg;
    }

    cmdline_parser_free (&args_info);
    double t[splinelength], x[splinelength], c_xy[splinelength];
    double * qrs_tmpl = malloc (splinelength*sizeof (double));
    unsigned int m_corr[splinelength/(2*qrs_length)];// = malloc (sizeof (unsigned int) * splinelength/200);
    unsigned int m_corr1[splinelength/(2*qrs_length)];
    unsigned int m_corr2[splinelength/(2*qrs_length)];
    int interval = 100;
    int _dat[splinelength];
    double sigma[qrs_length];
    const size_t p = 4;
    gsl_spline * spline;
    const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver * s;
    gsl_matrix * covar = gsl_matrix_alloc (p, p);
    double result[p][splinelength/(4*qrs_length)];
    double error[p + 2][splinelength/(2*qrs_length)];

    nfft_plan pl[p];
    for (i = 0; i < splinelength/400; i++) /* half length because of symmetry */
    {
        result[0][i] = 0.0;
        result[1][i] = 0.0;
        result[2][i] = 0.0;
        result[3][i] = 0.0;
    }
    for (n = channel_start; n < channel_start + channel_count; n++)
    {
        //switch (n)
        //{
        //    case 5:
        //        continue;
        //    case 134:
        //        continue;
        //    case 135:
        //        continue;
        //    case 136:
        //        continue;
        //}
        read_signals (_dat, 4, n, splinelength, filename);
        for (i = 0; i < splinelength; i++) 
        {
            x[i] = (double) _dat[i];
            t[i] = ((double)i)/samplingrate;
            if (i >= qrs_start && i < qrs_start + qrs_length)
                qrs_tmpl[i - qrs_start] =_dat[i];
            else if (i >= qrs_start)
                qrs_tmpl[i - qrs_start] = 0.0;
        }

        int scale = cspl_norm (x, splinelength);
        for (int l = 0; l < 2; l++)
        {
            cspl_norm (qrs_tmpl, qrs_length);
            a = 0;
            b = 0;
            k = 0;
            //if (args_info.enum_opt_given)
            {
                //  if (strcmp (args_info.enum_opt_arg, cmdline_parser_enum_opt_values[0]) == 0)
                {
                    cspl_xcorr (c_xy, qrs_tmpl, x, qrs_length, splinelength);
                    cspl_norm (c_xy, splinelength);
                    /*
                       for (i = 0; i < splinelength; i++) {
                       printf("%d %f\n", i, c_xy[i]);
                       }
                       */
                    a = cspl_eval_periodic_max2 (m_corr1, c_xy, splinelength, 2*qrs_length, 0.15);
                }
            }
            //else
            {
                cspl_root_mean_square (c_xy, qrs_tmpl, x, qrs_length, splinelength);
                cspl_norm (c_xy, splinelength);
                /* 
                   for (i = 0; i < splinelength; i++) {
                   printf("%d %f\n", i, c_xy[i]);
                   }
                   */
                b = cspl_eval_periodic_min2 (m_corr2, c_xy, splinelength, 2*qrs_length, -0.15);
            }
            //m_corr = malloc (sizeof (unsigned int) * GSL_MIN (a, b));
            for (i = 0; i < a; i++) {
                for (j = 0; j < b; j++)
                {
                    if ((m_corr1[i] == m_corr2[j]) || (m_corr1[i] == m_corr2[j] - 1) || (m_corr1[i] == m_corr2[j] + 1))
                    {
                        m_corr[k] = m_corr1[i];
                        k++;
                    }
                }
            }
            //m_corr[k] = 0;
            //printf ("#m=0,S=2\n");
            //for (i = 0; i < GSL_MAX(a,b); i++) 
            //{
            //    printf ("%d, %d\n", m_corr1[i], m_corr2[i]);
            //    printf ("%d %d\n",i, m_corr[i]);
            //}
            if (k == 0) 
            {
                continue;
            }
            {
                //qrs_tmpl = realloc (qrs_tmpl, splinelength/200*sizeof (double));
                for (i = 0; i < splinelength/(2*qrs_length); i++)
                {
                    qrs_tmpl[i] = 0.0;
                }
            }
            interval = cspl_average (qrs_tmpl, x, m_corr, splinelength, k);
        }
        if (interval > 700 || interval < 300){
            continue;
        }
        //interval = qrs_length; /* restrict the interval to the QRS-Komplex */

        //printf("n=%d, k=%d, interval=%d\n", n, k, interval);
        s = gsl_multifit_fdfsolver_alloc (T, qrs_length, p);
        //sigma = (double *) realloc (sigma, interval*sizeof (double));

        cspl_norm (qrs_tmpl, qrs_length);
        if (a > b)
        {
            k = a;
            memcpy (m_corr, m_corr1, k*sizeof (unsigned int));
        }
        else
        {
            k = b;
            memcpy (m_corr, m_corr2, k*sizeof (unsigned int));
        }
        if (k < 100)
        {
            continue;
        }
        if (k < k_min)
        {
            k_min = k;
        }
        //printf("interval=%d,k_min=%d\n", interval, k_min);
        cspl_qrs_init();
        spline = gsl_spline_alloc (gsl_interp_cspline, qrs_length);
        gsl_spline_init (spline, t, qrs_tmpl, qrs_length);
        if (template) {
            printf ("#m=1,S=0\n");
            for (j = 0; j < qrs_length*10; j++) /* spline with interpolation of 10 points */
                printf ("%.5f %.5f\n", t[j]/10, gsl_spline_eval(spline, t[j]/10, acc));
            //for (j = 0; j < interval; j++)
            //    printf ("%.5f %.5f\n", t[j], qrs_tmpl[j]);
        }

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
        double c0, c1, cov00, cov01, cov11, sumsq;
        double chisq_pdof, c;
        for (i = 0; i < k; i++) {
            double x_init[4] = { 1.0, 0.0, 0.0, 0.0 };
            gsl_vector_view vec = gsl_vector_view_array (x_init, p);
            gsl_fit_linear (t, 1, x + m_corr[i] + 400, 1, 150, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
            for (int m = 0; m < qrs_length; m++) {
                sigma[m] = sqrt(sumsq/150);
            }
            struct cspl_qrs_data data = {qrs_length, t, x + m_corr[i], spline, sigma, p, s, covar, vec, c, chisq_pdof};
            cspl_qrs_fit (&data);
            size_t j;
            pl[0].f[i] = gsl_vector_get(s->x, 0);
            pl[1].f[i] = gsl_vector_get(s->x, 1);
            pl[2].f[i] = gsl_vector_get(s->x, 2);
            pl[3].f[i] = gsl_vector_get(s->x, 3);
            error[0][i] = data.c*sqrt(gsl_matrix_get(covar,0,0));
            error[1][i] = data.c*sqrt(gsl_matrix_get(covar,1,1));
            error[2][i] = data.c*sqrt(gsl_matrix_get(covar,2,2));
            error[3][i] = data.c*sqrt(gsl_matrix_get(covar,3,3));
            error[4][i] = data.chisq_pdof;
            error[5][i] = sigma[0]*scale;
            if (i < 1)
            {
                if (timedomain)
                    printf ("#m=1,S=0\n");
                for (j = 0; j < qrs_length; j++) 
                {
                    double y_j = gsl_vector_get(s->x, 0)*gsl_spline_eval(spline, t[j] - gsl_vector_get(s->x, 1), acc) + gsl_vector_get(s->x, 2) + gsl_vector_get(s->x, 3)*(t[j] - gsl_vector_get(s->x, 1));
                    if (timedomain)
                        printf ("%.5f %.5f\n", t[j + m_corr[i]], y_j);
                }
                if (timedomain)
                    printf ("#m=2,S=0\n");
                for (j = m_corr[i]; j < m_corr[i] + qrs_length; j++) 
                {
                    if (timedomain)
                        printf ("%.5f %.5f\n", t[j], x[j]);
                }
            }
        }
        if (errors)
        {
            printf ("#m=2,S=0\n");
            for (i = 0; i < k; i++)
                printf("%.5f %.5f\n", t[m_corr[i]], error[4][i]);
            printf ("#m=2,S=0\n");

            for (i = 0; i < k; i++)
                printf("%.5f %.5f\n", t[m_corr[i]], error[5][i]);
        }
        //printf("mean(sigma)=%.5f, sd(sigma)=%.5f\n", gsl_stats_mean(error[5], 1, k), gsl_stats_sd(error[5], 1, k));
        gsl_multifit_fdfsolver_free (s);
        gsl_spline_free(spline);
        cspl_qrs_free ();


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

        //result = (double *) realloc (result, k*sizeof (double));
        for (i = 0; i < k/2; i++)
        {
            result[0][i] += cabs(pl[0].f_hat[i + k/2]);
            result[1][i] += cabs(pl[1].f_hat[i + k/2]);
            result[2][i] += cabs(pl[2].f_hat[i + k/2]);
            result[3][i] += cabs(pl[3].f_hat[i + k/2]);
        }
        /** finalise the one dimensional plan */
        nfft_finalize(&pl[0]);
        nfft_finalize(&pl[1]);
        nfft_finalize(&pl[2]);
        nfft_finalize(&pl[3]);
    }
    if (frq & 0x8)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k_min/2; i++) 
        {
            printf("%.5f %.5f\n", i/t[splinelength - 1], result[0][i]);
        }
    }
    if (frq & 0x4)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k_min/2; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1], result[1][i]);
        }
    }
    if (frq & 0x2)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k_min/2; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1], result[2][i]);
        }
    }
    if (frq & 0x1)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k_min/2; i++) {
            printf("%.5f %.5f\n", i/t[splinelength - 1], result[3][i]);
        }
    }

    if (time & 0x8)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f %.8f\n", t[m_corr[i]], creal(pl[0].f[i]), error[0][i]);
        }
    }
    if (time & 0x4)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f %.8f\n", t[m_corr[i]], creal(pl[1].f[i]), error[1][i]);
        }
    }
    if (time & 0x2)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f %.8f\n", t[m_corr[i]], creal(pl[2].f[i]), error[2][i]);
        }
    }
    if (time & 0x1)
    {
        printf ("#m=1,S=0\n");
        for (i = 0; i < k; i++) {
            printf("%.5f %.5f %.8f\n", t[m_corr[i]], creal(pl[3].f[i]), error[3][i]);
        }
    }
    gsl_matrix_free (covar);

    free (qrs_tmpl);
    free (filename);
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

