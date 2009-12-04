#include "cspl/cspl.h"

int signal_processing_eval_periodic_max (size_t * max_n, double * signal, size_t size, size_t period) {
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


int signal_processing_norm_average (double * signal, double * template, size_t * n, size_t size, size_t period) {
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
