/*   cspl_math.c                                                           *
 ***************************************************************************
 *   Copyright (C) 2009, 2010 Hermann Sonntag                              *
 *   hermann.sonntag@tu-ilmenau.de                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/
#include "cspl/cspl.h"

int cspl_root_mean_square (double * rms, double * x, double *y, size_t length, size_t size) {
    size_t i, j;
    for (i = 0; i < size; i++) {
        rms[i] = 0;
        if (i < size - length)
        for (j = i; j < i + length; j++) {
            rms[i] += pow(x[j - i] - y[j] - (x[0] - y[i]), 2);
        }
        else
        for (j = 0; j < length; j++) {
            rms[i] += pow(x[j] - y[j] - (x[0] - y[0]), 2);
        }

        if (length == 0)
            return GSL_EZERODIV;
        rms[i] /= length;
    }
    return GSL_SUCCESS;
}

int cspl_eval_periodic_min2 (unsigned int * min_n, double * signal, size_t size, size_t length, double min) {
    size_t i;
    size_t n = 0;
    double local_min = DBL_MAX;
    min_n[0] = n;
    for (i = 0; i < size - length - 1; i++) {
        while (signal[min_n[i]] > min) {
            for (n = min_n[i]; n < (min_n[i] + length) && (n < size); n++) {
                if (signal[n] < min) {
                    if (signal[n] < local_min) {
                        local_min = signal[n];
                        min_n[i] = n;
                        min_n[i + 1] = n + length;
                    }
                }
                else if (signal[min_n[i]] > min && n == min_n[i] + length - 1) {
                    min_n[i] = n;
                }
            }
        }
        if (signal[n] < min) {
            min_n[i] = n;
            min_n[i + 1] = n + length;
        }
        local_min = DBL_MAX;
        if(n == size)
            break;
    }
    return i;
}

int cspl_eval_periodic_min (unsigned int * min_n, double * signal, size_t size, double level) {
    size_t i;
    size_t n;
    int in_local_min = 0;
    double global_min = DBL_MAX;
    double local_min = DBL_MAX;
    min_n[0] = 0;

    for (i = 0; i < size - 1; i++) {
        for (n = min_n[i]; n < size; n++) {
            if ( signal[n] < global_min) 
                global_min = signal[n];
            if ( signal[n] < local_min) {
                local_min = signal[n];
                min_n[i] = n;
            }

            if ( signal[n] < global_min/level )
                in_local_min = 1;
            else if (in_local_min) {
                in_local_min = 0;
                local_min = DBL_MAX;
                min_n[i + 1] = n;
                break;
            }
        }
        if(n == size)
            break;
    }
    return i;
}
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
return i;
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
    return z;
}

int cspl_norm (double * signal, size_t size) {
    double max = DBL_MIN;
    double min = DBL_MAX;
    double mean = 0;
    size_t i;
    for (i = 0; i < size; i++) {
        mean += signal[i];
        if (signal[i] > max)
            max = signal[i];
        else if (signal[i] < min)
            min = signal[i];
    }
    if (size == 0)
        return GSL_EZERODIV;
    mean /= size;
    for (i = 0; i < size; i++) {
        signal[i] -= mean;
        if (max - min == 0)
            return GSL_EZERODIV;
        signal[i] /= (max - min);
    }
    return GSL_SUCCESS;
}

int cspl_average (double * templ, double * signal, unsigned int * n, size_t size, int count) {
    /*double max = DBL_MIN;
      double min = DBL_MAX;
      */
    unsigned int interval = UINT_MAX;
    size_t i, m;
    if (count == 0)
        return GSL_FAILURE;
    for (i = 0; i < count; i++) {
        if (n[i + 1] - n[i] < interval)
            interval = n[i + 1] - n[i];
    }
    for (i = 0; i < count; i++) {
        for (m = n[i]; m < n[i] + interval; m++) {
            templ[m - n[i]] += signal[m];
        }
    }
    /*for (i = 0; i < size; i++) {
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
      if (max - min == 0)
      return GSL_FAILURE;
      templ[i] /= (max - min);
      }
      */
    /*for (i = 0; i < size; i++) {
      if (i < interval) {
      templ[i] /= count;
      }
      else
      templ[i] = 0.0;
      }
      */
    return interval;
}
