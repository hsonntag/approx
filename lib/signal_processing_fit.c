#include "signal_processing/signal_processing.h"


gsl_interp_accel *acc;

struct signal_processing_QRS_data {
  size_t n;
  double * t;
  double * y;
  gsl_spline * n_QRS;
  double * sigma;
};
int signal_processing_QRS_init() {
acc = gsl_interp_accel_alloc();
return GSL_SUCCESS;
}
int signal_processing_QRS_free() {
  gsl_interp_accel_free(acc);
  return GSL_SUCCESS;
}

int signal_processing_radix2_xcorr (double * c_xy, double * x, double * y, size_t size) {
 int z;
  z = gsl_fft_real_radix2_transform(x, 1, size);
  z = gsl_fft_real_radix2_transform(y, 1, size);

  c_xy[0] = x[0]*y[0];
  c_xy[size/2] = x[size/2]*y[size/2];
  size_t i;
  for (i = 1; i < size/2; i++) {
    c_xy[i] = x[i]*y[i] + x[size-i]*y[size-i];
    c_xy[size-i] = x[size - i]*y[i] - x[i]*y[size-i];
  }
  z = gsl_fft_halfcomplex_radix2_inverse(c_xy, 1, size);
  return z;
}

int signal_processing_eval_max_periodic (size_t * max_n, double * signal, size_t size, size_t period) {
	size_t i;
	size_t n;
	double max;

	for (i = 0; i < size; i++) {
		max = 0;
		for (n = max_n[i]; (n < (max_n[i] + period)) && (n < size); n++) {
			if (signal[n] > max) {
				max = signal[n];
				max_n[i] = n;
			}
		}
		if((n == size) || (i == (size - 2)))
			break;
		max_n[i + 1] = max_n[i] +  (int)(0.5*period);
	}
return i;
}

int signal_processing_norm (double * signal, double * template, size_t * n, size_t size, size_t period) {
double max = 0;
size_t i, m;
for (i = 0; i < size; i++) {
	for (m = n[i]; m < n[i] + period; m++) {
		template[m - n[i]] += signal[m];
		if (template[m - n[i]]>max)
			max = template[m - n[i]];
	}
}

for (i = 0; i < size; i++) {
	template[i] /= max;
}
return GSL_SUCCESS;
}
/*------------------------*/
int signal_processing_fit (size_t n, double * t, double * y, gsl_spline * n_QRS, double * sigma) {
const gsl_multifit_fdfsolver_type * T;
gsl_multifit_fdfsolver * s;
int status;
unsigned int iter = 0;
//const size_t n = N;
const size_t p = 3;

gsl_matrix * covar = gsl_matrix_alloc (p, p);
struct signal_processing_QRS_data d = { n, t, y, n_QRS, sigma};
gsl_multifit_function_fdf f;
double x_init[3] = { 2.0, 0.0, 2.0 };
gsl_vector_view x = gsl_vector_view_array (x_init, p);
const gsl_rng_type * type;
gsl_rng * r;

gsl_rng_env_setup();

type = gsl_rng_default;
r = gsl_rng_alloc (type);

f.f = &signal_processing_QRS_f;
f.df = &signal_processing_QRS_df;
f.fdf = &signal_processing_QRS_fdf;
f.n = n;
f.p = p;
f.params = &d;

/* This is the data to be fitted */
size_t i;
for (i = 0; i < n; i++)
{
	/*  double t = i; */
	/*       y[i] = 1.0 + 5 * exp (-0.1 * t)  */
	/* 	+ gsl_ran_gaussian (r, 0.1); */
	sigma[i] = 0.1;
	/*  printf ("data: %u %g %g\n", i, y[i], sigma[i]); */
};

T = gsl_multifit_fdfsolver_lmsder;
s = gsl_multifit_fdfsolver_alloc (T, n, p);
gsl_multifit_fdfsolver_set (s, &f, &x.vector);

print_state (iter, s);

do
{
	iter++;
	status = gsl_multifit_fdfsolver_iterate (s);

	printf ("status = %s\n", gsl_strerror (status));

	print_state (iter, s);

	if (status)
		break;

	status = gsl_multifit_test_delta (s->dx, s->x,
			1e-4, 1e-4);
}
while (status == GSL_CONTINUE && iter < 500);

gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

{ 
	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

	printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof); 

	printf ("a      = %.5f +/- %.5f\n", FIT(0), c*ERR(0)); 
	printf ("t_beat = %.5f +/- %.5f\n", FIT(1), c*ERR(1)); 
	printf ("S_0      = %.5f +/- %.5f\n", FIT(2), c*ERR(2)); 
}

printf ("status = %s\n", gsl_strerror (status));

gsl_multifit_fdfsolver_free (s);
gsl_matrix_free (covar);
gsl_rng_free (r);
return GSL_SUCCESS;
}

int signal_processing_QRS_f (const gsl_vector * x, void * data, 
		gsl_vector * f) {

	size_t n = ((struct signal_processing_QRS_data *)data)->n;
	double *y = ((struct signal_processing_QRS_data *)data)->y; 
	double *n_QRS = ((struct signal_processing_QRS_data *)data)->n_QRS;
	double *t = ((struct signal_processing_QRS_data *)data)->t;
	double *sigma = ((struct signal_processing_QRS_data *) data)->sigma;

	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);

	gsl_spline_init(spline, t, n_QRS, n);

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);
	double s_0 = gsl_vector_get (x, 2);

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Model yi = a * n_QRS(t[i] - t_beat) + s_0 */

		double y_i = a * gsl_spline_eval(spline, t[i] - t_beat, acc) + s_0;
		gsl_vector_set (f, i, (y_i - y[i])/sigma[i]);
	}
	gsl_spline_free(spline);
	return GSL_SUCCESS;
}


int signal_processing_QRS_df (const gsl_vector * x, void * data, 
		gsl_matrix * J) {

	double a = gsl_vector_get (x, 0);
	double t_beat = gsl_vector_get (x, 1);

	size_t n = ((struct signal_processing_QRS_data *)data)->n;
	double *n_QRS = ((struct signal_processing_QRS_data *)data)->n_QRS;
	double *t = ((struct signal_processing_QRS_data *)data)->t;
	double *sigma = ((struct signal_processing_QRS_data *) data)->sigma;

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);

	gsl_spline_init(spline, t, n_QRS, n);

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (yi - yi_fit)/sigma[i],      */
		/*       yi = a * n_QRS(t[i] - t_beat) + S_0  */
		/* and the xj are the parameters (a,t_beat,S_0) */
//		((struct signal_processing_QRS_data *)data)->i = i;

		double s = sigma[i];
		double dQRS_da = gsl_spline_eval (spline, t[i] - t_beat, acc);
		double dQRS_dt_beat = -a*gsl_spline_eval_deriv (spline,  t[i] - t_beat,  acc);
		gsl_matrix_set (J, i, 0, dQRS_da/s); 
		gsl_matrix_set (J, i, 1, dQRS_dt_beat/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return GSL_SUCCESS;
}

int signal_processing_QRS_fdf (const gsl_vector * x, void *data,
		gsl_vector * f, gsl_matrix * J)
{
	signal_processing_QRS_f (x, data, f);
	signal_processing_QRS_df (x, data, J);

	return GSL_SUCCESS;
}


