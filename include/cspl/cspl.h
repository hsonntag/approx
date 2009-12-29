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

int cspl_radix2_xcorr(double * c_xy, double * x, double * y, size_t size);
int cspl_eval_periodic_max (unsigned int * max_n, double * signal, size_t size, double level);
int cspl_norm_average (double * templ, double * signal, unsigned int * n, size_t size, int count);

__END_DECLS

#endif /* __CSPL_H__ */
