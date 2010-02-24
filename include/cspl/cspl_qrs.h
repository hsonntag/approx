/*   cspl/cspl_qrs.h                                                       *
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
#ifndef __CSPL_QRS_H__
#define __CSPL_QRS_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlin.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

gsl_interp_accel * acc;

struct cspl_qrs_data {
    size_t n;
    double * t;
    double * y;
    gsl_spline * n_qrs;
    double * sigma;
    size_t p;
    gsl_multifit_fdfsolver * s;
    gsl_matrix * covar;
    gsl_vector_view x;
};

int cspl_qrs_init();

int cspl_qrs_free();

int cspl_qrs_fdf (const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);
int cspl_qrs_f (const gsl_vector * x, void * data, gsl_vector * f);
int cspl_qrs_df (const gsl_vector * x, void * data, gsl_matrix * J);

struct cspl_qrs_function_fdf_struct
{
    int (* f) (const gsl_vector * x, void * params, gsl_vector * f);
    int (* df) (const gsl_vector * x, void * params, gsl_matrix * df);
    int (* fdf) (const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix *df);
    size_t n;   /* number of functions */
    size_t p;   /* number of independent variables */
    void * params;
};

typedef struct cspl_qrs_function_fdf_struct cspl_qrs_function_fdf ;

__END_DECLS

#endif /* __CSPL_QRS_H__ */
