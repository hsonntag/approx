#include <signal_processing/signal_processing.h>
gsl_interp_accel *acc;

struct signal_processing_QRS_data {
  size_t n;
  double * t;
  double * y;
  gsl_spline * n_QRS;
  double * sigma;
};
int signal_processing_QRS_fdf (const gsl_vector * x, void * data,
		gsl_vector * f, gsl_matrix * J);
int signal_processing_QRS_f (const gsl_vector * x, void * data, gsl_vector * f);
int signal_processing_QRS_df (const gsl_vector * x, void * data, gsl_matrix * J);
