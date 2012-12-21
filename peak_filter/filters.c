#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "filters.h"
#include "numeric_constants.h"

double max(double *v, const unsigned int n)
{
    double max = DBL_MIN;
    unsigned int i;
    
    for (i = 0; i < n; i++)
    	if (v[i] > max) max = v[i];
    	
    return max;
}

unsigned int get_max_index(double *v, const unsigned int n)
{
    unsigned int i, max_index = 0;

    for (i = 1; i < n; i++)
    	if (v[i] > v[max_index]) max_index = i;

    return max_index;
}

double min(double *v, const unsigned int n)
{
    double min = DBL_MAX;
    unsigned int i;

    for (i = 0; i < n; i++)
    	if (v[i] < min) min = v[i];
    	
    return min;
}

unsigned int threshold(double *v, const unsigned int n, double threshold_value)
{
    unsigned int i, valid_elements = 0;

    for (i = 0; i < n; i++) {
    	if (v[i] >= threshold_value) {
    		v[valid_elements] = v[i];
    		valid_elements++;
        }
	}

    return valid_elements;
}


unsigned int threshold2(double *v, double *w, const unsigned int n, double threshold_value)
{
    unsigned int i, valid_elements = 0;

    for (i = 0; i < n; i++) {
    	if (v[i] >= threshold_value) {
    		v[valid_elements] = v[i];
      		w[valid_elements] = w[i];
    		valid_elements++;
        }
	}

    return valid_elements;
}

double sum(double *v, const unsigned int n)
{
    unsigned int i;
    double total = 0;

    for (i = 0; i < n; i++)
    	total += v[i];

    return total;
}

double average(double *v, const unsigned int n)
{
    return (sum(v, n) / n);
}

double cubic_average(double *v, const unsigned int n)
{
    unsigned int i;
    double total = 0;

    for (i = 0; i < n; i++)
    	total += pow(v[i], 0.33333333333333333333333333333333);

    return pow(total/n, 3);
}

double average2(double *v, const unsigned int n)
{
    double total = 0;
    unsigned int i, last_not_null_index;
    
    for (i = n - 1; v[i] == 0; i--);
    last_not_null_index = i;

    for (i = 0; i <= last_not_null_index; i++)
    	total += v[i];

    return total / (last_not_null_index + 1);
}

double average3(double *v, const unsigned int n)
{
    double total = 0;
    unsigned int i, first_not_null_index = 0, last_not_null_index;

    for (i = 0; v[i] == 0; i++);
    first_not_null_index = i;

    for (i = n - 1; v[i] == 0; i--);
    last_not_null_index = i;

    for (i = first_not_null_index; i <= last_not_null_index; i++)
    	total += v[i];

    return total / (last_not_null_index - first_not_null_index + 1);
}

double variance(double *v, const unsigned int n)
{
    double total_variance = 0, average_val = average(v, n);
    unsigned int i;

    for (i = 0; i < n; i++)
    	total_variance += pow((v[i] - average_val), 2);

    return (total_variance / n);
}

double variance2(double *v, const unsigned int n)
{
    double total_variance = 0, average_val = average(v, n);
    unsigned int i, total_i = 0;

    for (i = 0; i < n; i++) {
    	total_variance += (pow((v[i] - average_val), 2) * i);
    	total_i += i;
    }

    return (total_variance / total_i);
}


double weighted_variance(double *v, double *w, const unsigned int n)
{
    double total_variance = 0, average_val = average(v, n);
    unsigned int i;

    for (i = 0; i < n; i++)
    	total_variance += ( pow((v[i] - average_val), 2) * w[i] );

    return (total_variance / n);
}

int round_to_int(const double x)
{
	if (x >= 0) {
	    return (int)(x + 0.5);
	} else {
		return (int)(x - 0.5);
	}
}

unsigned int remove_masked(double amp[], double freq[], const unsigned int n,
											const double low_ratio, const double high_ratio)
{
    double *or_amp, *or_freq;
    unsigned int i, num_of_peaks;
    
    or_amp = (double*)malloc(n * sizeof(double));
    or_freq = (double*)malloc(n * sizeof(double));
    memcpy(or_amp, amp, n * sizeof(double));
    memcpy(or_freq, freq, n * sizeof(double));

    num_of_peaks = 0;
    while (1)
    {
    	unsigned int masker_index;
     	double masker_freq, masker_amp;
    
        masker_index = get_max_index(or_amp, n);
        masker_freq = or_freq[masker_index];
        masker_amp = or_amp[masker_index];
    
        if (masker_amp == DBL_MIN) break;
        
        freq[num_of_peaks] = masker_freq;
        amp[num_of_peaks] = masker_amp;
        num_of_peaks++;
        or_amp[masker_index] = DBL_MIN;
    
        for (i = 0; i < n; i++) {
            if (or_amp[i] != DBL_MIN
            	&& or_freq[i] >= masker_freq / low_ratio
            	&& or_freq[i] <= high_ratio * masker_freq)
             {
             	if (or_freq[i] <= masker_freq) {
             		if (or_amp[i] <= masker_amp * ( log(or_freq[i] * low_ratio / masker_freq)
               							/ log(low_ratio) ) )
                    {
                       or_amp[i] = DBL_MIN;
               		}
             	}
             	else {
                  	if (or_amp[i] <= masker_amp * ( (-1) * log(or_freq[i] / (high_ratio * masker_freq))
                  						/ log(high_ratio) ) )
                  	{
                       or_amp[i] = DBL_MIN;
                  	}
             	}
             };
        }
    }
    
    free(or_amp);
    free(or_freq);

    /* ordena por frequencia */
    for (i = num_of_peaks; i > 0; i--)
    {
    	unsigned int max_freq_index;
    	double aux;
    	
    	max_freq_index = get_max_index(freq, i);
    	
        aux = freq[i - 1];
        freq[i - 1] = freq[max_freq_index];
        freq[max_freq_index] = aux;
         
        aux = amp[i - 1];
        amp[i - 1] = amp[max_freq_index];
        amp[max_freq_index] = aux;
    }
    
    return num_of_peaks;
}


/*****************************************************************************************/
/* MASK_FILTER                                                                           */
/*****************************************************************************************/
unsigned int mask_filter(double amp[], double freq[], const unsigned int n,
											const double low_ratio, const double high_ratio)
{
    double *or_amp, *or_freq;
    unsigned int i, num_of_peaks;
    
    or_amp = (double*)malloc(n * sizeof(double));
    or_freq = (double*)malloc(n * sizeof(double));
    memcpy(or_amp, amp, n * sizeof(double));
    memcpy(or_freq, freq, n * sizeof(double));

    num_of_peaks = 0;
    while (1)
    {
    	unsigned int masker_index;
     	double masker_freq, masker_amp;
    
        masker_index = get_max_index(or_amp, n);
        masker_freq = or_freq[masker_index];
        masker_amp = or_amp[masker_index];
    
        if (masker_amp == DBL_MIN) break;
        
        freq[num_of_peaks] = masker_freq;
        amp[num_of_peaks] = masker_amp;
        num_of_peaks++;
        or_amp[masker_index] = DBL_MIN;
    
        for (i = 0; i < n; i++) {
            if (or_amp[i] != DBL_MIN
            	&& or_freq[i] >= masker_freq / low_ratio
            	&& or_freq[i] <= high_ratio * masker_freq)
             {
             	if (or_freq[i] <= masker_freq) {
             		if (or_amp[i] <= masker_amp * ( log(or_freq[i] * low_ratio / masker_freq)
               							/ log(low_ratio) ) )
                    {
                       or_amp[i] = DBL_MIN;
               		}
             	}
             	else {
                  	if (or_amp[i] <= masker_amp * ( (-1) * log(or_freq[i] / (high_ratio * masker_freq))
                  						/ log(high_ratio) ) )
                  	{
                       or_amp[i] = DBL_MIN;
                  	}
             	}
             };
        }
    }
    
    free(or_amp);
    free(or_freq);

    /* ordena por frequencia */
    for (i = num_of_peaks; i > 0; i--)
    {
    	unsigned int max_freq_index;
    	double aux;
    	
    	max_freq_index = get_max_index(freq, i);
    	
        aux = freq[i - 1];
        freq[i - 1] = freq[max_freq_index];
        freq[max_freq_index] = aux;
         
        aux = amp[i - 1];
        amp[i - 1] = amp[max_freq_index];
        amp[max_freq_index] = aux;
    }
    
    return num_of_peaks;
}

unsigned int remove_spurious_peaks(double amp[], double freq[], const unsigned int n)
{
    double *or_amp, *or_freq;
    unsigned int i, num_of_peaks;

    or_amp = (double*)malloc(n * sizeof(double));
    or_freq = (double*)malloc(n * sizeof(double));
    memcpy(or_amp, amp, n * sizeof(double));
    memcpy(or_freq, freq, n * sizeof(double));

    num_of_peaks = 0;
    while (1)
    {
    	unsigned int masker_index;
     	double masker_freq, masker_amp;

        masker_index = get_max_index(or_amp, n);
        masker_freq = or_freq[masker_index];
        masker_amp = or_amp[masker_index];

        if (masker_amp == DBL_MIN) break;

        freq[num_of_peaks] = masker_freq;
        amp[num_of_peaks] = masker_amp;
        num_of_peaks++;
        or_amp[masker_index] = DBL_MIN;

        for (i = 0; i < n; i++) {
            if (or_amp[i] != DBL_MIN)
            {
                if(freq[i] <= masker_freq)
                {
                    if (or_amp[i] <= (masker_amp - 60 - 18 * masker_freq / or_freq[i]) )
                    {
                        or_amp[i] = DBL_MIN;
                    }
                }
                else
                {
                    if (or_amp[i] <= (masker_amp - 60 - 18 * or_freq[i] / masker_freq) )
                    {
                        or_amp[i] = DBL_MIN;
                    }
                }
            }
        }
    }

    free(or_amp);
    free(or_freq);

    /* ordena por frequencia */
    for (i = num_of_peaks; i > 0; i--)
    {
    	unsigned int max_freq_index;
    	double aux;

    	max_freq_index = get_max_index(freq, i);

        aux = freq[i - 1];
        freq[i - 1] = freq[max_freq_index];
        freq[max_freq_index] = aux;

        aux = amp[i - 1];
        amp[i - 1] = amp[max_freq_index];
        amp[max_freq_index] = aux;
    }

    return num_of_peaks;
}

unsigned int new_filter_old(double *amp, double *freq, const unsigned int n,
																double max_per_octave_decay)
/* assumes non-decreasing freq order */
/* assumes amplitude in decibels */
{
    unsigned int i, j, max_index = get_max_index(amp, n);

    for (i = 0; i <= max_index; i++)
    {
        if (amp[i] <= 0) continue;
        
    	const double
     	 big_amp = amp[i],
      	 log_big_freq = log10(freq[i]);
    
        for (j = (i+1); j < n; j++)
        {
        	if (amp[j] <= 0) continue;

            const double actual_dif = big_amp - amp[j],
            		max_allowed_dif = max_per_octave_decay
              							* (log10(freq[j]) - log_big_freq)/LOG_OF_2_IN_BASE_10;
            
            if (actual_dif > max_allowed_dif) {
            	amp[j] = -1;
            }
        }
    }
    return threshold2(amp, freq, n, 0);
}

unsigned int new_filter(double *amp, double *freq, const unsigned int n,
																double max_per_harm_decay)
/* assumes non-decreasing freq order */
/* assumes amplitude in decibels */
{
    unsigned int i, j, max_index = get_max_index(amp, n);

    for (i = 0; i <= max_index; i++)
    {
        if (amp[i] <= 0) continue;

    	const double
     	 big_amp = amp[i],
      	 big_freq = freq[i];

        for (j = (i+1); j < n; j++)
        {
        	if (amp[j] <= 0) continue;

            const double actual_dif = big_amp - amp[j],
            		max_allowed_dif = max_per_harm_decay * (freq[j] / big_freq);

            if (actual_dif > max_allowed_dif) {
            	amp[j] = -1;
            }
        }
    }
    return threshold2(amp, freq, n, 0);
}


unsigned int amg_filter(double *amp, double *freq, const unsigned int n,
																double max_adj_harm_gap)
/* assumes non-decreasing freq order */
/* assumes amplitude in decibels */
{
    const unsigned int max_index = get_max_index(amp, n);
    unsigned int ref;
    signed int curr;
    
    for (ref = max_index, curr = (ref-1); curr >= 0; curr--)
    {
        double curr_ref_gap = amp[ref] - amp[curr];

    	if (curr_ref_gap > max_adj_harm_gap) {
        	amp[curr] = -1;
      	} else {
      		ref = curr;
      	}
    }
    
    for (ref = max_index, curr = (ref+1); curr < n; curr++)
    {
        double curr_ref_gap = amp[ref] - amp[curr];

    	if (curr_ref_gap > max_adj_harm_gap) {
        	amp[curr] = -1;
      	} else {
      		ref = curr;
      	}
    }
    
    return threshold2(amp, freq, n, 0);
}


/********************************************************************************************/
/* MULTI_FILTER                                                                             */
/********************************************************************************************/
/* does not implement the average filter */
unsigned int multi_filter(double *mag, double *freq, const unsigned int num_of_peaks,
							const double min_freq, const double max_freq,
							const double min_mag,
							const double max_mag_gap,
							const double max_adj_mag_gap)
{
	unsigned int i, strongest_index, curr_num_of_peaks;
	double strongest_mag, threshold_value;

	if (num_of_peaks == 0) {
		return 0;
	}
	
	for (i = 0, strongest_mag = mag[0], strongest_index = 0;
		 i < num_of_peaks;
		 i++)
	{
		if (freq[i] < min_freq || mag[i] < min_mag || freq[i] > max_freq) {
			mag[i] = DBL_MIN;
		} else {
			if (mag[i] > strongest_mag) {
				strongest_mag = mag[i];
				strongest_index = i;
			}
		}
		
	}
	
	if (strongest_mag < min_mag) {
		return 0;
	}

	{
		unsigned int valid_elements;
		threshold_value = strongest_mag - max_mag_gap;
		
	    for (i = 0, valid_elements = 0; i < num_of_peaks; i++) {
	    	if (mag[i] == DBL_MIN) continue;
	    	if (mag[i] >= threshold_value) {
	    		mag[valid_elements] = mag[i];
	      		freq[valid_elements] = freq[i];

	      		if (i == strongest_index) {
	      			strongest_index = valid_elements;
	      		}

	    		valid_elements++;
	        }
		}
		
		curr_num_of_peaks = valid_elements;
	}
	
	if (curr_num_of_peaks <= 1) {
		return curr_num_of_peaks;
	}

	{
		signed int ref, curr;
	    
	    for (ref = strongest_index, curr = (ref-1); curr >= 0; curr--)
	    {
	        double curr_ref_gap = mag[ref] - mag[curr];
	
	    	if (curr_ref_gap > max_adj_mag_gap) {
	        	mag[curr] = DBL_MIN;
	      	} else {
	      		ref = curr;
	      	}
	    }
	    
	    for (ref = strongest_index, curr = (ref+1); curr < curr_num_of_peaks; curr++)
	    {
	        double curr_ref_gap = mag[ref] - mag[curr];
	
	    	if (curr_ref_gap > max_adj_mag_gap) {
	        	mag[curr] = DBL_MIN;
	      	} else {
	      		ref = curr;
	      	}
	    }
	}

	{
		unsigned int valid_elements;
		
	    for (i = 0, valid_elements = 0; i < curr_num_of_peaks; i++) {
	    	if (mag[i] == DBL_MIN) continue;
	    	if (mag[i] >= threshold_value) {
	    		mag[valid_elements] = mag[i];
	      		freq[valid_elements] = freq[i];

	    		valid_elements++;
	        }
		}
		return valid_elements;
	}
}
