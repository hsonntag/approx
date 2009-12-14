#include "cspl/cspl.h"

int cspl_eval_periodic_max (size_t * max_n, double * signal, size_t size, size_t period) {
	size_t i;
	size_t n;
	double max;
	max_n[0] = 0;

	for (i = 0; i < size; i++) {
		max = 0;
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
	double tmp;
/*	for (i = 1; i < size/2; i++) {
		tmp = c_xy[i];
		c_xy[i] = c_xy[size - i];
		c_xy[size - i] = tmp;
	}
*/	return z;
}

int cspl_norm_average (double * templ, double * signal, size_t * n, size_t size, int count) {
	double max = -10000000.0;
	double min =  10000000.0;
	size_t i, m;
	for (i = 0; i < count; i++) {
		if (i == size - 2)
			break;
		for (m = n[i]; m < n[i + 1]; m++) {
			templ[m - n[i]] += signal[m];
			if (templ[m - n[i]] > max)
				max = templ[m - n[i]];
			else if (templ[m -n[i]] < min)
				min = templ[m - n[i]];
		}
	}

	for (i = 0; i < size; i++) {
		templ[i] -= min;
		templ[i] /= (max - min);
	}
	return GSL_SUCCESS;
}
