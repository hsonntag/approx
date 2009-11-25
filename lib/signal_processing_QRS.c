#include "signal_processing/signal_processing_QRS.h"

int signal_processing_QRS_init() {
acc = gsl_interp_accel_alloc();
return GSL_SUCCESS;
}
int signal_processing_QRS_free() {
  gsl_interp_accel_free(acc);
  return GSL_SUCCESS;
}

int signal_processing_QRS_f (const gsl_vector * x, void * data, 
		gsl_vector * f) {

	size_t n = ((struct signal_processing_QRS_data *)data)->n;
	double * y = ((struct signal_processing_QRS_data *)data)->y; 
	gsl_spline * n_QRS = ((struct signal_processing_QRS_data *)data)->n_QRS;
	double * t = ((struct signal_processing_QRS_data *)data)->t;
	double * sigma = ((struct signal_processing_QRS_data *) data)->sigma;


//	gsl_spline_init(spline, t, n_QRS, n);

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);
	double s_0 = gsl_vector_get (x, 2);

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Model yi = a * n_QRS(t[i] - t_beat) + s_0 */

		double y_i = a * gsl_spline_eval(n_QRS, t[i] - t_beat, acc) + s_0;
		gsl_vector_set (f, i, (y_i - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}


int signal_processing_QRS_df (const gsl_vector * x, void * data, 
		gsl_matrix * J) {

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);

	size_t n = ((struct signal_processing_QRS_data *)data)->n;
	gsl_spline * n_QRS = ((struct signal_processing_QRS_data *)data)->n_QRS;
	double * t = ((struct signal_processing_QRS_data *)data)->t;
	double * sigma = ((struct signal_processing_QRS_data *) data)->sigma;

	//gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, n);

//	gsl_spline_init(spline, t, n_QRS, n);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double s = sigma[i];
		double dQRS_da = gsl_spline_eval (n_QRS, t[i] - t_beat, acc);
		double dQRS_dt_beat = -a*gsl_spline_eval_deriv (n_QRS,  t[i] - t_beat,  acc);
		gsl_matrix_set (J, i, 0, dQRS_da/s); 
		gsl_matrix_set (J, i, 1, dQRS_dt_beat/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	return GSL_SUCCESS;
}

int signal_processing_QRS_fdf (const gsl_vector * x, void *data,
		gsl_vector * f, gsl_matrix * J)
{
	signal_processing_QRS_f (x, data, f);
	signal_processing_QRS_df (x, data, J);

	return GSL_SUCCESS;
}
