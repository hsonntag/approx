#include "cspl/cspl.h"

int cspl_eval_periodic_max (unsigned int * max_n, double * signal, size_t size, double level) {
	size_t i;
	size_t n;
	int in_local_max = 0;
	double global_max = DBL_MIN;
	double local_max = DBL_MIN;
	max_n[0] = 0;

	for (i = 0; i < size - 1; i++) {
		for (n = max_n[i]; n < size; n++) {
			if ( signal[n] > global_max) 
				global_max = signal[n];
			if ( signal[n] > local_max) {
				local_max = signal[n];
				max_n[i] = n;
			}

			if ( signal[n] > level * global_max )
				in_local_max = 1;
			else if (in_local_max) {
				in_local_max = 0;
				local_max = DBL_MIN;
				max_n[i + 1] = n;
				break;
			}
		}
		if(n == size)
			break;
	}
	/*
	   for (n = max_n[i]; (n < (max_n[i] + period/2)) && (n < size); n++) {
	   if (signal[n] > max) {
	   printf("i=%d, signal[%d] = %f\n", i, n, signal[n]);
	   max = signal[n];
	   max_n[i] = n;
	   }
	   }
	   if((n == size) || (i == (size - 2)))
	   break;
	   max_n[i + 1] = max_n[i] + period/2;
	   }
	   */
return i + 1;
}

int cspl_radix2_xcorr (double * c_xy, double * x, double * y, size_t size) {
	size_t i;
	int z;
	double _x[size], _y[size];
	for (i = 0; i < size; i++) {
		_x[i] = x[i];
		_y[i] = y[i];
	}
	z = gsl_fft_real_radix2_transform(_x, 1, size);
	z = gsl_fft_real_radix2_transform(_y, 1, size);

	c_xy[0] = _x[0]*_y[0];
	c_xy[size/2] = _x[size/2]*_y[size/2];
	for (i = 1; i < size/2; i++) {
		c_xy[i] = _x[i]*_y[i] + _x[size-i]*_y[size-i];
		c_xy[size - i] = _x[size - i]*_y[i] -_x[i]*_y[size-i];
	}
	z = gsl_fft_halfcomplex_radix2_inverse(c_xy, 1, size);
/*	double tmp;
		for (i = 1; i < size/2; i++) {
		tmp = c_xy[i];
		c_xy[i] = c_xy[size - i];
		c_xy[size - i] = tmp;
		}
*/		return z;
}

int cspl_norm (double * signal, size_t size) {
	double max = DBL_MIN;
	double min = DBL_MAX;
	size_t i;
	for (i = 0; i < size; i++) {
		if (signal[i] > max)
			max = signal[i];
		else if (signal[i] < min)
			min = signal[i];
	}
	for (i = 0; i < size; i++) {
		signal[i] -= min;
		signal[i] /= (max - min);
	}
	return GSL_SUCCESS;
}

int cspl_norm_average (double * templ, double * signal, unsigned int * n, size_t size, int count) {
	double max = DBL_MIN;
	double min = DBL_MAX;
	unsigned int interval = UINT_MAX;
	size_t i, m;
	for (i = 0; i < count; i++) {
		if (i == size - 2)
			break;

		if (n[i + 1] - n[i] < interval)
			interval = n[i + 1] - n[i];

		for (m = n[i]; m < n[i + 1]; m++) {
			templ[m - n[i]] += signal[m];
		}
	}

	for (i = 0; i < size; i++) {
		if (i < interval) {
			if (templ[i] > max)
				max = templ[i];
			else if (templ[i] < min)
				min = templ[i];
		}
		else 
			templ[i] = 0.0; 
	}
	for (i = 0; i < interval; i++) {
		templ[i] -= min;
		templ[i] /= (max - min);
	}
	return interval;
}
