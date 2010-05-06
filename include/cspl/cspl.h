/*   cspl/cspl.h                                                           *
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
#ifndef __CSPL_H__
#define __CSPL_H__

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
#include <cspl/cspl_qrs.h>
#include <cspl/cspl_qrs_fit.h>

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

int cspl_real_fft (double * data, size_t size);
int cspl_halfcomplex_abs (double * data, size_t size);
int cspl_root_mean_square (double * rms, double * x, double * y, size_t length, size_t size);
int cspl_norm (double * signal, size_t size);
int cspl_eval_periodic_min2 (unsigned int * min_n, double * signal, size_t size, size_t length, double min); 
int cspl_eval_periodic_max2 (unsigned int * max_n, double * signal, size_t size, size_t length, double max); 
int cspl_radix2_xcorr(double * c_xy, double * x, double * y, size_t size);
int cspl_xcorr(double * c_xy, double * x, double * y, size_t length, size_t size);
int cspl_eval_periodic_max (unsigned int * max_n, double * signal, size_t size, double level);
int cspl_average (double * templ, double * signal, unsigned int * n, size_t size, int count);

__END_DECLS

#endif /* __CSPL_H__ */
