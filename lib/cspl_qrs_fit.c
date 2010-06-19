/*   cspl_qrs_fit.c                                                        *
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
#include "cspl/cspl_qrs.h"
#include "cspl/cspl_qrs_fit.h"

int cspl_qrs_fit (void * params) {
    int status;
    unsigned int iter;
    struct cspl_qrs_data * data = (struct cspl_qrs_data *) params;
    /* This is the data to be fitted */
    size_t i;
    for (i = 0; i < data->n; i++)
    {
        data->sigma[i] = 1.0;
    }

    gsl_multifit_function_fdf f;
    //    const gsl_rng_type * type;
    //    gsl_rng * r;

//    gsl_rng_env_setup();

    //    type = gsl_rng_default;
    //    r = gsl_rng_alloc (type);

    f.f = &cspl_qrs_f;
    f.df = &cspl_qrs_df;
    f.fdf = &cspl_qrs_fdf;
    f.n = data->n;
    f.p = data->p;
    f.params = data;

    for (int j = 0; j < 2; j++)
    {

        gsl_multifit_fdfsolver_set (data->s, &f, &data->x.vector);
        iter = 0;
        //print_state (iter, data->s);

        do
        {
            iter++;
            status = gsl_multifit_fdfsolver_iterate (data->s);
#ifdef DEBUG
            printf ("status = %s\n", gsl_strerror (status));

            print_state (iter, data->s);
#endif
            if (status)
                break;

            status = gsl_multifit_test_delta (data->s->dx, data->s->x,
                    1e-12, 1e-12);
        }
        while (status == GSL_CONTINUE && iter < 500);

        gsl_multifit_covar (data->s->J, 0.0, data->covar);



        double chi = gsl_blas_dnrm2(data->s->f);
        double dof = data->n - data->p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
        double sigma = 0;
        for (i = 400; i < data->n; i++)
        {
            sigma += pow(gsl_vector_get(data->s->f, i), 2.0);
        }
        sigma = sqrt(sigma/(dof - 400));
        if (j == 0)
        {
            for (i = 0; i < data->n; i++)
            {
                if (sigma != 0.0)
                    data->sigma[i] = sigma;
            }
        }
        else
        {
            data->sigma[0] = pow(chi, 2.0)/dof;
            data->sigma[1] = c;
        }
    }
#ifdef DEBUG
#define FIT(i) gsl_vector_get(data->s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(data->covar,i,i))
    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof); 

    printf ("A        = %.5f +/- %.5f\n", FIT(0), c*ERR(0)); 
    printf ("t_beat   = %.5f +/- %.5f\n", FIT(1), c*ERR(1)); 
    printf ("S_0      = %.5f +/- %.5f\n", FIT(2), c*ERR(2)); 
    printf ("S_1      = %.5f +/- %.5f\n", FIT(3), c*ERR(3)); 


    printf ("status = %s\n", gsl_strerror (status));
#endif
    //    gsl_rng_free (r);
    return GSL_SUCCESS;
}

