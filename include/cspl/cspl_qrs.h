#ifndef __CSPL_QRS_H__
#define __CSPL_QRS_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>

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

gsl_interp_accel *acc;

struct cspl_qrs_data {
	size_t n;
	double * t;
	double * y;
	gsl_spline * n_qrs;
	double * sigma;
};

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
