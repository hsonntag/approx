#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>


struct data {
  size_t n;
  size_t i;
  double * t;
  double * y;
  double * n_QRS;
  double * sigma;
};
     
int
s_QRS_f (const gsl_vector * x, void *data, 
	 gsl_vector * f)
{

  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y; 
  double *n_QRS = ((struct data *)data)->n_QRS;
  double *t = ((struct data *)data)->t;
  double *sigma = ((struct data *) data)->sigma;
     
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);

  gsl_spline_init(spline, t, n_QRS, n);

  double A = gsl_vector_get (x, 0);
  double t_beat = gsl_vector_get (x, 1);
  double S_0 = gsl_vector_get (x, 2);
     
  size_t i;
     
  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * n_QRS(t[i] - t_beat) + S_0 */

      double Yi = A * gsl_spline_eval(spline, t[i] - t_beat, acc) + S_0;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return GSL_SUCCESS;
}

/* double */
/* n_QRS_t_beat (double t_beat, void *data) */
/* { */

/*   size_t n = ((struct data *)data)->n; */
/*   size_t i = ((struct data *)data)->i; */
/*   double *n_QRS = ((struct data *)data)->n_QRS; */
/*   double *t = ((struct data *)data)->t; */
     
/*   gsl_interp_accel *acc = gsl_interp_accel_alloc(); */
/*   gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n); */

/*   gsl_spline_init(spline, t, n_QRS, n); */

/*   double Yi = gsl_spline_eval(spline, t[i] - t_beat, acc); */
 
/*   gsl_spline_free(spline); */
/*   gsl_interp_accel_free(acc); */
/*   return Yi; */
/* }  */

/* double */
/* n_QRS_t_i (double t_i, void *data) */
/* { */

/*   size_t n = ((struct data *)data)->n; */
/*   double *n_QRS = ((struct data *)data)->n_QRS; */
/*   double *t = ((struct data *)data)->t; */
     
/*   gsl_interp_accel *acc = gsl_interp_accel_alloc(); */
/*   gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n); */

/*   gsl_spline_init(spline, t, n_QRS, n); */

/*   double Yi = gsl_spline_eval(spline, t_i, acc); */
 
/*   gsl_spline_free(spline); */
/*   gsl_interp_accel_free(acc); */
/*   return Yi; */
/* }   */ 
 
int
s_QRS_df (const gsl_vector * x, void *data, 
	  gsl_matrix * J)
{
     
  double A = gsl_vector_get (x, 0);
  double t_beat = gsl_vector_get (x, 1);
     
  // gsl_function F;
  //double result, abserr;
     
  //F.function = &n_QRS_t_beat;
  //F.params = data;

  size_t n = ((struct data *)data)->n;
  // double *y = ((struct data *)data)->y; 
  double *n_QRS = ((struct data *)data)->n_QRS;
  double *t = ((struct data *)data)->t;
  double *sigma = ((struct data *) data)->sigma;
     
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);

  gsl_spline_init(spline, t, n_QRS, n);

  size_t i;
     
  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * n_QRS(t[i] - t_beat) + S_0  */
      /* and the xj are the parameters (A,t_beat,S_0) */
      ((struct data *)data)->i = i;
      //gsl_deriv_central (&F, t_beat, 1e-8, &result, &abserr);

      double s = sigma[i];
      double a = gsl_spline_eval (spline, t[i] - t_beat, acc);
      double result = -gsl_spline_eval_deriv (spline,  t[i] - t_beat,  acc);
      gsl_matrix_set (J, i, 0, a/s); 
      gsl_matrix_set (J, i, 1, A*result/s);
      gsl_matrix_set (J, i, 2, 1/s);
    }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return GSL_SUCCESS;
}
     
int
s_QRS_fdf (const gsl_vector * x, void *data,
	   gsl_vector * f, gsl_matrix * J)
{
  s_QRS_f (x, data, f);
  s_QRS_df (x, data, J);
     
  return GSL_SUCCESS;
}


