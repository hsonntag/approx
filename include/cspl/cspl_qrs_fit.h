/*   cspl/cspl_qrs_fit.h                                                   *
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
#ifndef __CSPL_QRS_FIT_H__
#define __CSPL_QRS_FIT_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
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

//extern int print_state (unsigned int iter, char * s);
int cspl_qrs_fit (double * t, double * y, gsl_spline * n_qrs, double * sigma, size_t n);
int cspl_qrs_fit_at (double * t, double * y, gsl_spline * n_qrs, double * sigma, size_t n, int count);

__END_DECLS

#endif /* __CSPL_QRS_FIT_H__ */
