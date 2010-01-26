#include "cspl/cspl_qrs.h"

int cspl_qrs_init() {
	acc = gsl_interp_accel_alloc();
	return GSL_SUCCESS;
}

int cspl_qrs_free() {
	gsl_interp_accel_free(acc);
	return GSL_SUCCESS;
}

int cspl_qrs_f (const gsl_vector * x, void * data, gsl_vector * f) {

	size_t n = ((struct cspl_qrs_data *)data)->n;
	double * y = ((struct cspl_qrs_data *)data)->y; 
	gsl_spline * n_qrs = ((struct cspl_qrs_data *)data)->n_qrs;
	double * t = ((struct cspl_qrs_data *)data)->t;
	double * sigma = ((struct cspl_qrs_data *) data)->sigma;


	//	gsl_spline_init(spline, t, n_qrs, n);

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);
	double s_0 = gsl_vector_get (x, 2);

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Model yi = a * n_qrs(t[i] - t_beat) + s_0 */

		double y_i = a * gsl_spline_eval(n_qrs, t[i] - t_beat, acc) + s_0;
		printf("f(%d) = %f * N_QRS(%f - %f) + %f = %f\n", i, a, t[i], t_beat, s_0, y_i);
		gsl_vector_set (f, i, (y_i - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}


int cspl_qrs_df (const gsl_vector * x, void * data, 
		gsl_matrix * J) {

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);

	size_t n = ((struct cspl_qrs_data *)data)->n;
	gsl_spline * n_qrs = ((struct cspl_qrs_data *)data)->n_qrs;
	double * t = ((struct cspl_qrs_data *)data)->t;
	double * sigma = ((struct cspl_qrs_data *) data)->sigma;

	//gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, n);

	//	gsl_spline_init(spline, t, n_qrs, n);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double s = sigma[i];
		double dqrs_da = gsl_spline_eval (n_qrs, t[i] - t_beat, acc);
		double dqrs_dt_beat = -a*gsl_spline_eval_deriv (n_qrs,  t[i] - t_beat,  acc);
		gsl_matrix_set (J, i, 0, dqrs_da/s); 
		gsl_matrix_set (J, i, 1, dqrs_dt_beat/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	return GSL_SUCCESS;
}

int cspl_qrs_fdf (const gsl_vector * x, void *data,
		gsl_vector * f, gsl_matrix * J)
{
	cspl_qrs_f (x, data, f);
	cspl_qrs_df (x, data, J);

	return GSL_SUCCESS;
}

cspl_qrs_function_fdf cspl_qrs_fnc = {
	.f = cspl_qrs_f,
	.df = cspl_qrs_df,
	.fdf = cspl_qrs_fdf,
};
