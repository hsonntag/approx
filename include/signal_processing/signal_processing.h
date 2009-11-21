#include <string.h>
#include <math.h>
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
//double * x, * y;
//size_t size;
int signal_processing_radix2_xcorr(double * c_xy, double * x, double * y, size_t size);
int signal_processing_QRS_fdf (const gsl_vector * x, void * data,
		gsl_vector * f, gsl_matrix * J);
int signal_processing_QRS_f (const gsl_vector * x, void * data, gsl_vector * f);
int signal_processing_QRS_df (const gsl_vector * x, void * data, gsl_matrix * J);
