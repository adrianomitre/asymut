#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "f0_estimation_methods.h"
#include "data_structures.h"
#include "numeric_constants.h"
#include "config_constants.h"
#include "aux_fun.h"


static double *critical_band_weight = NULL;


double average_f0_est(double freq[], double amp[], const unsigned int n)
{
    unsigned int i;
    double total_freq = 0, total_amp = 0;

    for (i = 0; i < n; i++)
    {
    	total_freq += freq[i] / (i+1) * amp[i];
    	total_amp += amp[i];
    }

    return total_freq / total_amp;
}

double cb_average_f0_est(double freq[], double amp[], const unsigned int n,
																	const double cb_weight)
/* new f0 estimation criteria, which uses the critical bands concept */
{
    unsigned int i;
    double total_freq = 0, total_amp = 0, weight;

    for (i = 0; i < n; i++)
    {
    	weight = 1.0 / pow((i + 1), cb_weight);
    	total_freq += freq[i] / (i+1) * amp[i] * weight;
    	total_amp += amp[i] * weight;
    }

    return total_freq / total_amp;
}


/********************************************************************************************/
/* FULL_AVERAGE_F0_EST                                                                      */
/********************************************************************************************/
double full_average_f0_est(double freq[], double amp[], const unsigned int n,
	const double min_f0, const double max_f0,
    	const double cb_weight_error, const double cb_weight_f0_est)
/* expects NON DECREASING freq[], or will return -1 */
{
    const double total_energy = sum(amp, n);
    
    unsigned int i, j, max_harm, best_f0_index, num_of_candidates, num_of_accum_cand;
    double *candidates /* TODO: refactor as cand_freq */, **cand_freq, **cand_amp,
    		possible_candidate, highest_freq_peak, strongest_peak_freq,
      		curr_f0_est, curr_f0_error, best_f0_est, min_f0_error,
        	max_acceptable_f0_candidate = DBL_MAX, min_candidate_energy, accum_energy;
        	
    if (n == 1) {
        if (freq[0] >= min_f0 && freq[0] <= max_f0) {
            return freq[0];
        }
        else {
        	return 0;
        }
    }

    min_candidate_energy = sum(amp, n) * MIN_CANDIDATE_TO_TOTAL_RATIO;
    for (i = 0, accum_energy = 0; i < n; i++) {
        accum_energy += amp[i];
        if (accum_energy >= MAX_PERCENTILE_FREQ_CANDIDATE * total_energy) {
        	max_acceptable_f0_candidate = (1.0 / MIN_HARM_TO_FUND_RATIO) * freq[i];
            break;
        }
    }

    candidates = (double*)emalloc((n + 1) * sizeof(double));
    cand_freq = (double**)emalloc((n + 1) * sizeof(double*));
    cand_amp = (double**)emalloc((n + 1) * sizeof(double*));
    
    highest_freq_peak = freq[n - 1];
    strongest_peak_freq = freq[max_index(amp, n)];
   {
    num_of_candidates = 0;
    num_of_accum_cand = 0;
    for (i = 0; i < (n - 1); i++) {
        possible_candidate = freq[i + 1] - freq[i];

        if (possible_candidate < 0) {
        	return -1;
        }
        
        if (possible_candidate <= 2 * strongest_peak_freq
        	&& possible_candidate <= max_acceptable_f0_candidate)
        {
        	candidates[num_of_candidates] = possible_candidate;
            num_of_candidates++;
        }

    }
    
    candidates[num_of_candidates++] = freq[0];
    candidates[num_of_candidates++] = strongest_peak_freq;
   }

    max_harm = (highest_freq_peak / min(candidates, num_of_candidates) + 0.5);

    for (i = 0; i < num_of_candidates; i++) {
        cand_freq[i] = (double*)emalloc(max_harm * sizeof(double));
        cand_amp[i] = (double*)emalloc(max_harm * sizeof(double));
        
        for (j = 0; j < max_harm; j++) {
        	cand_freq[i][j] = DBL_MAX;
        	cand_amp[i][j] = 0;
        }
    }
    
    for (i = 0; i < num_of_candidates; i++) {
        for (j = 0; j < n; j++) {
            if (freq[j] <= MIN_HARM_TO_FUND_RATIO * candidates[i]) continue;
            
            unsigned int closest_harm = freq[j] / candidates[i] - 0.5;
            /*
            double expected_freq = (closest_harm + 1) * candidates[i];
            */
            
            if (
            	/*
            	absoluteRelativeError(expected_freq, freq[j])
            	 <
           		absoluteRelativeError(expected_freq, cand_freq[i][closest_harm])
           		*/
           		amp[j] > cand_amp[i][closest_harm]
               )
            {
            	cand_freq[i][closest_harm] = freq[j];
                cand_amp[i][closest_harm] = amp[j];
            }
            
        }
    }

   {
    best_f0_est = 0;
    min_f0_error = DBL_MAX;
    best_f0_index = 0;

    for (i = 0; i < num_of_candidates; i++) {

    	if (max(cand_amp[i], max_harm) != max(amp, n)) continue;

        spec_smooth(cand_freq[i], cand_amp[i], max_harm);
        curr_f0_est = cb_average_f0_est(&cand_freq[i][0], &cand_amp[i][0], max_harm, cb_weight_f0_est);
        curr_f0_error = f0_cb_error(curr_f0_est, &cand_freq[i][0], &cand_amp[i][0], max_harm, cb_weight_error);
        /*
        curr_f0_est = average_f0_est(&cand_freq[i][0], &cand_amp[i][0], max_harm);
        curr_f0_error = f0_error(curr_f0_est, &cand_freq[i][0], &cand_amp[i][0], max_harm);
        */

        if (curr_f0_error < min_f0_error) {
        	if (curr_f0_est >= min_f0
                 &&
                curr_f0_est <= max_f0
        		 &&
        		(
        		min_f0_error == DBL_MAX
        		 ||
        		average2(cand_amp[i], max_harm) >= MIN_AVERAGE2_RATIO * average2(cand_amp[best_f0_index], max_harm)
                )
                &&
                sum(cand_amp[i], max_harm) >= min_candidate_energy
                &&
        		odd_average2(cand_amp[i], max_harm) >= MIN_ODD_AVERAGE2_RATIO * odd_average2(cand_amp[best_f0_index], max_harm)
           	   )
         	{
            	best_f0_index = i;
                best_f0_est = curr_f0_est;
                min_f0_error = curr_f0_error;
            }
        }
    }
   }

    free(candidates);
    for (i = 0; i < num_of_candidates; i++) {
        free(cand_freq[i]);
        free(cand_amp[i]);
    }
    free(cand_freq);
    free(cand_amp);

    return best_f0_est;
}


double f0_error(double f0, double freq[], double amp[], const unsigned int n)
{
	unsigned int i;
    double total_error = 0, expected_freq;

    for (i = 0; i < n; i++) {
        expected_freq = f0 * (i + 1);
        if (freq[i] <= expected_freq) {
        	total_error += (1 - freq[i]/expected_freq) * amp[i];
        }
        else {
        	total_error += (expected_freq/freq[i] - 1) * amp[i];
        }
    }

    return total_error;
}

/********************************************************************************************/
/* AVERAGE_ABS_F0_ERROR                                                                     */
/********************************************************************************************/
double averageAbsF0Error(double f0, double freq[], double amp[], const double totalEnergy, const unsigned int n)
{
	unsigned int i;
    double total_error = 0;
    signed int last_not_null_index = (n-1);
    
    while (amp[last_not_null_index] == 0 && last_not_null_index >= 0) {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 0; i <= last_not_null_index; i++)
    {
    	if (amp[i] > 0) {
	    	double expected_freq = f0 * (i + 1);

	        total_error += absoluteRelativeError(freq[i], expected_freq) * amp[i];
    	}
    }

    return total_error / (totalEnergy * (last_not_null_index + 1));
}


/********************************************************************************************/
/* TOTAL_ABS_F0_ERROR                                                                       */
/********************************************************************************************/
double totalAbsF0Error(double f0, double freq[], double amp[], const double totalEnergy, const unsigned int n)
{
	unsigned int i;
    double total_error = 0;
    signed int last_not_null_index = (n-1);
    
    while (amp[last_not_null_index] == 0 && last_not_null_index >= 0) {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 0; i <= last_not_null_index; i++)
    {
    	if (amp[i] > 0) {
	    	double expected_freq = f0 * (i + 1);

	        total_error += absoluteRelativeError(freq[i], expected_freq) * amp[i];
    	}
    }

    return total_error;
}


double f0_cb_error(double f0, double freq[], double amp[], const unsigned int n,
																	const double cb_weight)
/* new error measure, which uses the critical bands concept */
{
	unsigned int i;
    double total_error = 0, expected_freq, weigth;

    for (i = 0; i < n; i++) {
        if (freq[i] == 0 || amp[i] == 0) {
            continue;
        }

        weigth = 1.0 / pow((i + 1), cb_weight);
        expected_freq = f0 * (i + 1);
        if (freq[i] <= expected_freq) {
        	total_error += (expected_freq/freq[i] - 1) * amp[i] * weigth;
        }
        else {
        	total_error += (freq[i]/expected_freq - 1) * amp[i] * weigth;
        }
    }

    return total_error;
}


double max(double *v, const unsigned int n)
{
    double max;
    unsigned int i;
    
    for (i = 0, max = v[0]; i < n; i++)
    	if (v[i] > max) max = v[i];
    	
    return max;
}

unsigned int max_index(double *v, const unsigned int n)
{
    unsigned int i, max_index = 0;

    for (i = 1; i < n; i++)
    	if (v[i] > v[max_index]) max_index = i;

    return max_index;
}

double min(double *v, const unsigned int n)
{
    double min;
    unsigned int i;

    for (i = 0, min = v[0]; i < n; i++)
    	if (v[i] < min) min = v[i];
    	
    return min;
}

unsigned int min_index(double *v, const unsigned int n)
{
    unsigned int i, min_index = 0;

    for (i = 1; i < n; i++)
    	if (v[i] < v[min_index]) min_index = i;

    return min_index;
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
    double total = 0;
    unsigned int i;

    for (i = 0; i < n; i++)
    	total += v[i];

    return total;
}

/********************************************************************************************/
/* IS_PRIME                                                                                 */
/********************************************************************************************/
unsigned char is_prime(int n)
{
	unsigned int i;
	
	if (n < 0) {
		return is_prime(-n);
	}
	if (n <= 1) {
		return FALSE;
	}
	
	for (i = 2; i <= floor(sqrt(n) + 0.5); i++) {
		if (n%i == 0) {
			return FALSE;
		}
	}
	
	return TRUE;
}

/********************************************************************************************/
/* PREPARE_CRITICAL_BAND_SUM                                                                */
/********************************************************************************************/
void prepare_critical_band_sum(unsigned int n)
{
	unsigned int i;
	double h;
	
	critical_band_weight = (double*) emalloc(n * sizeof(double));

	critical_band_weight[0] = 1;
	for (i = 1, h = 2; i < n; i++, h += 1)
	{
		double a, b;

		a = sqrt(h/(h-1))*(h-1);
		b = sqrt((h+1)/h)*h;

		critical_band_weight[i] = minimum(1, log10(b)/THIRD_OCTAVE_UP_LOG10 - log10(a)/THIRD_OCTAVE_UP_LOG10);
		
		/*
		if (i > 1 && !is_prime(i)) {
			critical_band_weight[i] /= i;
		}
		*/
	}
}

/********************************************************************************************/
/* KLAP_PREPARE_CRITICAL_BAND_SUM                                                           */
/********************************************************************************************/
void klap_prepare_critical_band_sum(unsigned int n)
{
	unsigned int i;
	
	critical_band_weight = (double*) emalloc(n * sizeof(double));
	
	for (i = 0; i < n; i++)
	{
		double a, b;

		a = (i+1)*THIRD_OCTAVE_DOWN;
		b = (i+1)*THIRD_OCTAVE_UP;

		critical_band_weight[i] = 0.25 + 0.75/ceil(b-a+1);
	}
}

/********************************************************************************************/
/* MY_OLD_PREPARE_CRITICAL_BAND_SUM                                                         */
/********************************************************************************************/
void my_old_prepare_critical_band_sum(unsigned int n)
{
	unsigned int i;
	
	critical_band_weight = (double*) emalloc(n * sizeof(double));
	
	critical_band_weight[0] = 1;
	for (i = 2; i <= n; i++)
	{
		double a, b;

		a = sqrt((double)i/(i-1)) * (i-1);		
		b = sqrt((double)(i+1)/i) * i;

		critical_band_weight[i-1] = minimum(1, 1.0/(log10((b-a)/a)/THIRD_OCTAVE_LOG10));
	}
}

/********************************************************************************************/
/* FINALIZE_CRITICAL_BAND_SUM                                                               */
/********************************************************************************************/
void finalize_critical_band_sum(void)
{
	free(critical_band_weight);
}

/********************************************************************************************/
/* CRITICAL_BAND_SUM                                                                        */
/********************************************************************************************/
double critical_band_sum(double *v, const unsigned int n)
{
    unsigned int i;
    double total;

    for (i = 0, total = 0; i < n; i++) {
    	total += (critical_band_weight[i] * v[i]);
    }
    
    return total;
}

double penalize(const double absRelativeError)
/* expects non-negative values */
{
	/*
	const double ALPHA = 2.151/HALFTONE_UP, EXPONENT = 2;
	return exp(-pow(absRelativeError*ALPHA, EXPONENT));
	*/
	/*
	return 1-absRelativeError/HALFTONE_UP;
	*/
	return 1;
}

/********************************************************************************************/
/* NEW_CRITICAL_BAND_SUM                                                                    */
/********************************************************************************************/
double new_critical_band_sum(const double f0, double freq[], double mag[], const unsigned int n)
{
    unsigned int i;
    double total;

    for (i = 0, total = 0; i < n; i++) {
    	total += (critical_band_weight[i] * mag[i] * penalize(absoluteRelativeError(f0*i, freq[i])));
    }
    
    return total;
}

/********************************************************************************************/
/* ERB_SUM                                                                                  */
/********************************************************************************************/
double erb_sum(double freq[], double mag[], const unsigned int n)
{
	const unsigned char num_of_bands = 40; /* number of ERB bands */

    unsigned int k, b, band_num_of_bins[num_of_bands];
    double band_energy[num_of_bands], total_volume;

    for (b = 0; b < num_of_bands; b++)
    {
    	band_num_of_bins[b] = 0;
    	band_energy[b] = 0;
    }
    
    for (k = 0; k < n; k++)
    {
        b = floor(hertz_to_erb(freq[k]) + 0.5);
        
        if (b >= num_of_bands) break;

       	band_num_of_bins[b]++;
       	
        band_energy[b] += pow(10,mag[k]);
   	}
   	
	for (b = 0, total_volume = 0; b < num_of_bands; b++) {
        if (band_num_of_bins[b] >= 1) {
   	    	total_volume += log10((double)band_energy[b]/band_num_of_bins[b]);
      	}
	}

    return 10*total_volume; /* Bel to decibel */
}


/********************************************************************************************/
/* BARK_SUM                                                                                  */
/********************************************************************************************/
double bark_sum(double freq[], double mag[], const unsigned int n)
{
	const unsigned char num_of_bands = 26; /* number of Bark bands */

    unsigned int k, b, band_num_of_bins[num_of_bands];
    double band_energy[num_of_bands], total_volume;

    for (b = 0; b < num_of_bands; b++)
    {
    	band_num_of_bins[b] = 0;
    	band_energy[b] = 0;
    }
    
    for (k = 0; k < n; k++)
    {
        b = floor(hertz_to_bark(freq[k]) + 0.5);
        
        if (b >= num_of_bands) break;

       	band_num_of_bins[b]++;
       	
        band_energy[b] += pow(10,mag[k]);
   	}
   	
	for (b = 0, total_volume = 0; b < num_of_bands; b++) {
        if (band_num_of_bins[b] >= 1) {
   	    	total_volume += log10((double)band_energy[b]/band_num_of_bins[b]);
      	}
	}

    return 10*total_volume; /* Bel to decibel */
}


/********************************************************************************************/
/* CRITICAL_BAND_AVERAGE_F0_ESTIMATE                                                        */
/********************************************************************************************/
double critical_band_average_f0_estimate(double freq[], double amp[], const unsigned int n)
{
    unsigned int i;
    double total_freq, total_weight = 0;

    for (i = 0; i < n; i++)
    {
    	total_freq += freq[i] / (i+1) * amp[i] * critical_band_weight[i];
    	total_weight += amp[i] * critical_band_weight[i];
    }

    return total_freq / total_weight;
}


/********************************************************************************************/
/* CRITICAL_BAND_AVERAGE2                                                                   */
/********************************************************************************************/
double critical_band_average2(double *v, const unsigned int n)
{
    unsigned int i;
    double totalAmp, totalWeight;
    signed int last_not_null_index = (n-1);
    
    while (v[last_not_null_index] == 0 && last_not_null_index >= 0)	 {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 0, totalAmp = 0; i <= last_not_null_index; i++) {
    	totalAmp += (critical_band_weight[i] * v[i]);
    	totalWeight += critical_band_weight[i];
    }
    
    return totalAmp / totalWeight;
}


/********************************************************************************************/
/* UNSMOOTHNESS                                                                             */
/********************************************************************************************/
double unsmoothness(double *v, const unsigned int n)
{
    unsigned int i;
    double unsmoothness, totalWeight;
    signed int last_not_null_index = (n-1);
    
    while (v[last_not_null_index] == 0 && last_not_null_index >= 0)	 {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 1, unsmoothness = totalWeight = 0; i <= last_not_null_index; i++) {
    	double w = log((double)(i+1.0)/i);
    	totalWeight+= w;
    	unsmoothness += fabs(v[i]-v[i-1]) * w;
    }
    
    return unsmoothness/totalWeight;
}

double odd_even_energy_ratio(double *v, const unsigned int n)
{
}

/********************************************************************************************/
/* ODD_CRITICAL_BAND_AVERAGE2                                                               */
/********************************************************************************************/
double odd_critical_band_average2(double *v, const unsigned int n)
{
    unsigned int i;
    double totalAmp, totalWeight;
    signed int last_not_null_index = (n-1);
    
    while (v[last_not_null_index] == 0 && last_not_null_index >= 0)	 {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 0, totalAmp = 0; i <= last_not_null_index; i+=2) {
    	totalAmp += (critical_band_weight[i] * v[i]);
    	totalWeight += critical_band_weight[i];
    }
    
    return totalAmp / totalWeight;
}


double odd_sum(double *v, const unsigned int n)
{
    double total = 0;
    unsigned int i;

    for (i = 0; i < n; i += 2)
    	total += v[i];

    return total;
}

double average(double *v, const unsigned int n)
{
    return (sum(v, n) / n);
}

double average2(double *v, const unsigned int n)
{
    double total = 0;
    unsigned int i;
    signed int last_not_null_index = (n-1);
    
    while (v[last_not_null_index] == 0 && last_not_null_index >= 0)	 {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;
    
    for (i = 0; i <= last_not_null_index; i++) {
    	total += v[i];
    }

    return total / (last_not_null_index + 1);
}

double odd_average2(double *v, const unsigned int n)
{
    double total = 0;
    unsigned int i;
    signed int last_not_null_index = (n-1);
    
    while (v[last_not_null_index] == 0 && last_not_null_index >= 0)	 {
    	last_not_null_index--;
    }
    
    if (last_not_null_index < 0) return 0;

    for (i = 0; i <= last_not_null_index; i += 2) {
    	total += v[i];
    }

    return total / (last_not_null_index/2 + 1);
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
    return (int)(x + 0.5);
}

/********************************************************************************************/
/* SPEC_SMOOTH                                                                              */
/********************************************************************************************/
void spec_smooth(double freq[], double amp[], const unsigned int n)
{
    const double WIDTH = 2, LOG_WIDTH = log(WIDTH); /* in octaves */
    unsigned int i, k, last_index;
    double total_amp, total_weight, smoothed_amp,
    	curr_freq, curr_amp, min_freq, max_freq,
    	neighbour_freq, neighbour_amp, neighbour_weight;

    for (i = round_to_int(WIDTH); i < n; i++)
    {
    	curr_freq = freq[i];
    	curr_amp = amp[i];
    	
    	if (curr_freq <= 0 || curr_amp <= 0) {
     		continue;
        }
    	
        min_freq = curr_freq / WIDTH;
        max_freq = curr_freq * WIDTH;
        total_amp = 0;
        total_weight = 0;
        last_index = round_to_int((i+1.5) * WIDTH) - 1;
        if (last_index >= n) {
            last_index = n - 1;
        }

        for (k = round_to_int((i+0.5) / WIDTH) - 1; k <= last_index; k++)
        {
            neighbour_freq = freq[k];
            neighbour_amp = amp[k];

            if (neighbour_freq < min_freq || neighbour_freq > max_freq) {
                continue;
            }
            
            if (neighbour_freq <= curr_freq) {
            	neighbour_weight = log(neighbour_freq/min_freq)/LOG_WIDTH;
            }
            else {
                neighbour_weight = -log(neighbour_freq/max_freq)/LOG_WIDTH;
            }
            
            total_amp += neighbour_weight * neighbour_amp;
            total_weight += neighbour_weight;
        }
        
        smoothed_amp = total_amp / total_weight;
        
        if (smoothed_amp <= curr_amp) {
            /*
            fprintf(stderr, "amp: %g -> %g\n", amp[i], smoothed_amp);
            */
            amp[i] = smoothed_amp;
        }
    }

}


double absoluteRelativeError(const double f1, const double f2)
{
    if (f1 >= f2) {
        return (f1 - f2)/f2;
    } else {
        return (f2 - f1)/f1;
    }
}



double getFreqByAmpQuantileBelow(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks, const double quantile)
{
	const double quantAmp = quantile * sum(peaksAmp, numOfPeaks);
	unsigned int i;
	double accAmp;
	
	if (quantile == 1) return DBL_MAX;
	
	for (i = 0, accAmp = 0; i < numOfPeaks; i++) {
		accAmp += peaksAmp[i];
		if (accAmp > quantAmp) break;
	}
	
	return peaksFreq[i-1];
}

double getFreqByAmpQuantileAbove(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks, const double quantile)
{
	const double quantAmp = quantile * sum(peaksAmp, numOfPeaks);
	int i;
	double accAmp;
	
	if (quantile == 1) return 0;
	
	for (i = numOfPeaks-1, accAmp = 0; i > 0; i--) {
		accAmp += peaksAmp[i];
		if (accAmp > quantAmp) break;
	}
	
	return peaksFreq[i+1];
}

/********************************************************************************************/
/* MITRE_F0_ESTIMATE                                                                        */
/********************************************************************************************/
double mitre_f0_est(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks,
	const double time, const EstimationParameters *estPar)
/* expects NON DECREASING freq[], or will return -1 */
{
	unsigned int i, j, numOfCandidates, numOfValidCandidates, maxHarmNumber, F0EstIndex;
	unsigned char *isCandValid;
	double minCandFreq, maxCandFreq, F0EstFreq,
			*candGrossEst, *candGrossSalience, /* (antigo) usar a GrossSalience para obter F0 pela mediana */
			**candHarmFreq, **candHarmAmp,     /* (recente) REMOVER GrossSalience, pois parece ser inútil */
			*candAverageParcialEnergy, *candTotalEnergy, *candFreqError,
			strongestPeakFreq, minCandEnergy;
			
	if (estPar->opMode == PRINT_PARAMETERS_OM) {
		printf("OperationMode: %d\n", estPar->opMode);
		printf("minAbsFreq: %g\n", estPar->minAbsFreq);
		printf("maxAbsFreq: %g\n", estPar->maxAbsFreq);
		printf("maxRelativeEnergyBelow: %g\n", estPar->maxRelativeEnergyBelow);
		printf("minRelativeEnergy: %g\n", estPar->minRelativeEnergy);

		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* I. Produce initial constraints                                                        */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	/*
	minCandFreq = estPar->minAbsFreq;
	maxCandFreq = estPar->maxAbsFreq;
	*/
	minCandFreq = getFreqByAmpQuantileAbove(peaksFreq, peaksAmp, numOfPeaks, estPar->maxRelativeEnergyAbove);
	if (estPar->minAbsFreq > minCandFreq) {
		minCandFreq = estPar->minAbsFreq;
	}
	maxCandFreq	= getFreqByAmpQuantileBelow(peaksFreq, peaksAmp, numOfPeaks, estPar->maxRelativeEnergyBelow);
	if (estPar->maxAbsFreq < maxCandFreq) {
		maxCandFreq = estPar->maxAbsFreq;
	}
	if (estPar->opMode == PRINT_MAX_CAND_FREQ_OM) {
		printf("%g %g\n", time, maxCandFreq);
		
		return OK;
	}
	if (estPar->opMode == PRINT_MIN_CAND_FREQ_OM) {
		printf("%g %g\n", time, minCandFreq);
		
		return OK;
	}
	if (estPar->maxAbsFreq < maxCandFreq) maxCandFreq = estPar->maxAbsFreq;
	if (estPar->opMode == PRINT_MAX_CAND_FREQ_OM) {
		printf("%g %g\n", time, maxCandFreq);
		
		return OK;
	}

	if (maxCandFreq < minCandFreq) {
		return ERROR;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* II. Generate possible candidates                                                      */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	candGrossEst = (double*) emalloc(1024 * sizeof(double));
	candGrossSalience = (double*) emalloc(1024 * sizeof(double));
	
	strongestPeakFreq = peaksFreq[max_index(peaksAmp, numOfPeaks)];
	numOfCandidates = 0;

	for (i = 0; ;i++) {
		double possibleCand = strongestPeakFreq/(double)(i+1);

		if (possibleCand > maxCandFreq) {
			continue;
		}
		
		if (possibleCand >= minCandFreq) {
			candGrossEst[numOfCandidates] = possibleCand;
			candGrossSalience[numOfCandidates] = 1.0/(i+1);
			numOfCandidates++;
		} else {
			break;
		}
	}

	/*
	candGrossEst = (double*) emalloc((2*numOfPeaks - 1) * sizeof(double));
	candGrossSalience = (double*) emalloc((2*numOfPeaks - 1) * sizeof(double));

	for (i = 0, numOfCandidates = 0; i < numOfPeaks; i++)
	{
		if (peaksFreq[i] >= minCandFreq
			&& peaksFreq[i] <= maxCandFreq)
		{
			candGrossEst[numOfCandidates] = peaksFreq[i];
			candGrossSalience[numOfCandidates] = peaksAmp[i];
			numOfCandidates++;
		}
		
		if (i < (numOfPeaks - 1)) {
			double peaksFreqInterval = peaksFreq[i+1] - peaksFreq[i];
			
			if (peaksFreqInterval >= minCandFreq
				&& peaksFreqInterval <= maxCandFreq)
			{
				candGrossEst[numOfCandidates] = peaksFreqInterval;
				candGrossSalience[numOfCandidates] = 0.5*(peaksAmp[i+1] + peaksAmp[i]);
				numOfCandidates++;
			}
		}
	}
	*/

	if (numOfCandidates == 0) {
		return ERROR;
	}

	candGrossEst = (double*) erealloc(candGrossEst, numOfCandidates * sizeof(double));
	candGrossSalience = (double*) erealloc(candGrossSalience, numOfCandidates * sizeof(double));
	
	if (estPar->opMode == PRINT_INIT_CANDIDATES) {
		for (i = 0; i < numOfCandidates; i++) {
			printf("%g %g %g\n", time, candGrossEst[i], candGrossSalience[i]);
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* III. Group candidate partials (harmonic series)                                       */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	candHarmFreq = (double**) emalloc(numOfCandidates * sizeof(double*));
	candHarmAmp = (double**) emalloc(numOfCandidates * sizeof(double*));

	/* TODO: maxHarmNumber doubt: ceil or round? */
	maxHarmNumber = ceil(peaksFreq[numOfPeaks - 1] / min(candGrossEst, numOfCandidates));

	for (i = 0; i < numOfCandidates; i++) {
		candHarmFreq[i] = (double*) emalloc(maxHarmNumber * sizeof(double));
		candHarmAmp[i] = (double*) emalloc(maxHarmNumber * sizeof(double));
	}

	for (i = 0; i < numOfCandidates; i++) {
		for (j = 0; j < maxHarmNumber; j++) {
			candHarmFreq[i][j] = 0;
			candHarmAmp[i][j] = 0;
		}
	}

	candAverageParcialEnergy = (double*) emalloc(numOfCandidates * sizeof(double));
	candTotalEnergy = (double*) emalloc(numOfCandidates * sizeof(double));
	candFreqError = (double*) emalloc(numOfCandidates * sizeof(double));
	isCandValid = (unsigned char*) emalloc(numOfCandidates * sizeof(unsigned char));
	
	for (i = 0; i < numOfCandidates; i++) {
		for (j = 0; j < numOfPeaks; j++)
		{
			const double harm_ratio = peaksFreq[j]/candGrossEst[i];
			double low, high, inharm;
			unsigned int closest_harm;
			
			/* ROUND TO LOGARITHMICALLY CLOSEST HARMONIC */
			
			if (harm_ratio/floor(harm_ratio) <= ceil(harm_ratio)/harm_ratio) {
				closest_harm = floor(harm_ratio);
			} else {
				closest_harm = ceil(harm_ratio);
			}
			
			if (harm_ratio <= closest_harm) {

				/* removing the following 5 lines (if-else block) provokes bizarre errors! */
				if (closest_harm == 1) {
					low = 0;
				} else {
					low = (closest_harm-1.0)*sqrt(closest_harm/(closest_harm-1.0));
				}

				inharm = closest_harm/harm_ratio - 1;
			} else {
				/* removing the following line provokes bizarre errors! */
				high = closest_harm*sqrt((closest_harm+1.0)/closest_harm);
				
				inharm = harm_ratio/closest_harm - 1;
			}
			if (inharm <= estPar->maxRelativeInharmonicity
				&& (peaksAmp[j] > candHarmAmp[i][closest_harm-1]
					||
					(peaksAmp[j] == candHarmAmp[i][closest_harm-1]
						&& inharm < absoluteRelativeError(candHarmFreq[i][closest_harm-1], candGrossEst[i]*closest_harm))))
			{
				candHarmFreq[i][closest_harm-1] = peaksFreq[j];
				candHarmAmp[i][closest_harm-1] = peaksAmp[j];
			}
		}

		candTotalEnergy[i] = critical_band_sum(candHarmAmp[i], maxHarmNumber);
	}

	minCandEnergy = estPar->minRelativeEnergy * max(candTotalEnergy, numOfCandidates);

	if (estPar->opMode == PRINT_INIT_CAND_TOTAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			printf("%g %g %g\n", time, candGrossEst[i], candTotalEnergy[i]);
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* IV. Filter out invalid candidates                                                     */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	for (i = 0, numOfValidCandidates = numOfCandidates; i < numOfCandidates; i++) {
		if (candTotalEnergy[i] < minCandEnergy)	{
			isCandValid[i] = NO;

			candTotalEnergy[i] = 0;
			candAverageParcialEnergy[i] = DBL_MAX;
			
			numOfValidCandidates--;
		}
		else
		{
			isCandValid[i] = YES;
			/*
			candFreqError[i] = f0_error(average_f0_est(candHarmFreq[i], candHarmAmp[i], maxHarmNumber), candHarmFreq[i], candHarmAmp[i], maxHarmNumber);
			*/

			/*
			candAverageParcialEnergy[i] = average2(candHarmAmp[i], maxHarmNumber);
			*/
			/*
			spec_smooth(candHarmFreq[i], candHarmAmp[i], maxHarmNumber);
			*/
			/*
			candAverageParcialEnergy[i] = critical_band_average2(candHarmAmp[i], maxHarmNumber);
			*/
			candAverageParcialEnergy[i] = unsmoothness(candHarmAmp[i], maxHarmNumber);
			candFreqError[i] = averageAbsF0Error(candGrossEst[i], candHarmFreq[i], candHarmAmp[i], candTotalEnergy[i], maxHarmNumber);
		}
	}

	if (numOfValidCandidates == 0) {
		return ERROR;
	}
	
	if (estPar->opMode == PRINT_VALID_CAND_GROSS_FREQ_AND_TOTAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candTotalEnergy[i]);
			}
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	
	if (estPar->opMode == PRINT_VALID_CAND_GROSS_FREQ_AND_AVERAGE_PARTIAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candAverageParcialEnergy[i]);
			}
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	if (estPar->opMode == PRINT_VALID_CAND_ERROR) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candFreqError[i]);
			}
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* V. Rank candidates                                                                    */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	F0EstIndex = min_index(candAverageParcialEnergy, numOfCandidates);
	/*
	printf("candAverageParcialEnergy[F0EstIndex]=%g\n", candAverageParcialEnergy[F0EstIndex]);
	*/
	
	if (estPar->opMode == PRINT_F0_HARMONIC_SERIES) {
		for (j = 0; j < maxHarmNumber; j++) {
			if (candHarmAmp[F0EstIndex][j]) {
				printf("%g %g %g %u\n", time, candHarmFreq[F0EstIndex][j], candHarmAmp[F0EstIndex][j], j+1);
			}
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
			
	
	if (estPar->opMode == PRINT_GROSS_F0_ESTIMATE)
	{
		printf("%g %g %g\n", time, candGrossEst[F0EstIndex], candAverageParcialEnergy[F0EstIndex]);

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* VI. Refine frequency estimate                                                         */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	/*
	F0EstFreq = average_f0_est(candHarmFreq[F0EstIndex], candHarmAmp[F0EstIndex], maxHarmNumber);
	*/
	F0EstFreq = critical_band_average_f0_estimate(candHarmFreq[F0EstIndex], candHarmAmp[F0EstIndex], maxHarmNumber);


	/* TODO: usar spec_smooth() */

	if (estPar->opMode == PRINT_REFINED_F0_ESTIMATE)
	{
		if (F0EstFreq > 0) {
			/*
			fprintf(stderr, "DEBUG: %g %g --> %g\n", time, candGrossEst[F0EstIndex], F0EstFreq);
			*/
			printf("%g %g %g\n", time, F0EstFreq, candTotalEnergy[F0EstIndex]);
		} else {
			fprintf(stderr, "%g %g\n", time, F0EstFreq);
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* VII. Free allocated memory                                                            */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(candGrossEst);
	free(candGrossSalience);
	free(candAverageParcialEnergy);
	free(candTotalEnergy);
	free(candFreqError);
	free_2d((void**)candHarmFreq, numOfCandidates);
	free_2d((void**)candHarmAmp, numOfCandidates);


	return F0EstFreq;
}


/********************************************************************************************/
/* ERB_F0_ESTIMATE                                                                        */
/********************************************************************************************/
double erb_f0_est(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks,
	const double time, const EstimationParameters *estPar)
/* expects NON DECREASING freq[], or will return -1 */
{
	unsigned int i, j, numOfCandidates, numOfValidCandidates, maxHarmNumber, F0EstIndex;
	unsigned char *isCandValid;
	double minCandFreq, maxCandFreq, F0EstFreq,
			*candGrossEst, *candGrossSalience, /* (antigo) usar a GrossSalience para obter F0 pela mediana */
			**candHarmFreq, **candHarmAmp,     /* (recente) REMOVER GrossSalience, pois parece ser inútil */
			*candAverageParcialEnergy, *candTotalEnergy, *candFreqError,
			strongestPeakFreq, minCandEnergy;
	
	if (estPar->opMode == PRINT_PARAMETERS_OM) {
		printf("OperationMode: %d\n", estPar->opMode);
		printf("minAbsFreq: %g\n", estPar->minAbsFreq);
		printf("maxAbsFreq: %g\n", estPar->maxAbsFreq);
		printf("maxRelativeEnergyBelow: %g\n", estPar->maxRelativeEnergyBelow);
		printf("minRelativeEnergy: %g\n", estPar->minRelativeEnergy);

		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* I. Produce initial constraints                                                        */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	/*
	minCandFreq = estPar->minAbsFreq;
	maxCandFreq = estPar->maxAbsFreq;
	*/
	minCandFreq = getFreqByAmpQuantileAbove(peaksFreq, peaksAmp, numOfPeaks, estPar->maxRelativeEnergyAbove);
	if (estPar->minAbsFreq > minCandFreq) {
		minCandFreq = estPar->minAbsFreq;
	}
	maxCandFreq	= getFreqByAmpQuantileBelow(peaksFreq, peaksAmp, numOfPeaks, estPar->maxRelativeEnergyBelow);
	if (estPar->maxAbsFreq < maxCandFreq) {
		maxCandFreq = estPar->maxAbsFreq;
	}
	if (estPar->opMode == PRINT_MAX_CAND_FREQ_OM) {
		printf("%g %g\n", time, maxCandFreq);
		
		return OK;
	}
	if (estPar->opMode == PRINT_MIN_CAND_FREQ_OM) {
		printf("%g %g\n", time, minCandFreq);
		
		return OK;
	}
	if (estPar->maxAbsFreq < maxCandFreq) maxCandFreq = estPar->maxAbsFreq;
	if (estPar->opMode == PRINT_MAX_CAND_FREQ_OM) {
		printf("%g %g\n", time, maxCandFreq);
		
		return OK;
	}

	if (maxCandFreq < minCandFreq) {
		return ERROR;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* II. Generate possible candidates                                                      */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	candGrossEst = (double*) emalloc(MAX_F0_CANDIDATES * sizeof(double));
	candGrossSalience = (double*) emalloc(MAX_F0_CANDIDATES * sizeof(double));
	
	strongestPeakFreq = peaksFreq[max_index(peaksAmp, numOfPeaks)];
	
	for (i = 0, numOfCandidates = 0; numOfCandidates < MAX_F0_CANDIDATES;i++) {
		double possibleCand = strongestPeakFreq/(double)(i+1);

		if (possibleCand > maxCandFreq) {
			continue;
		}
		
		if (possibleCand >= minCandFreq) {
			candGrossEst[numOfCandidates] = possibleCand;
			candGrossSalience[numOfCandidates] = 1.0/(i+1);
			numOfCandidates++;
		} else {
			break;
		}
	}

	/*
	candGrossEst = (double*) emalloc((2*numOfPeaks - 1) * sizeof(double));
	candGrossSalience = (double*) emalloc((2*numOfPeaks - 1) * sizeof(double));

	for (i = 0, numOfCandidates = 0; i < numOfPeaks; i++)
	{
		if (peaksFreq[i] >= minCandFreq
			&& peaksFreq[i] <= maxCandFreq)
		{
			candGrossEst[numOfCandidates] = peaksFreq[i];
			candGrossSalience[numOfCandidates] = peaksAmp[i];
			numOfCandidates++;
		}
		
		if (i < (numOfPeaks - 1)) {
			double peaksFreqInterval = peaksFreq[i+1] - peaksFreq[i];
			
			if (peaksFreqInterval >= minCandFreq
				&& peaksFreqInterval <= maxCandFreq)
			{
				candGrossEst[numOfCandidates] = peaksFreqInterval;
				candGrossSalience[numOfCandidates] = 0.5*(peaksAmp[i+1] + peaksAmp[i]);
				numOfCandidates++;
			}
		}
	}
	*/

	if (numOfCandidates == 0) {
		return ERROR;
	}
	
	if (numOfCandidates != MAX_F0_CANDIDATES) {
		candGrossEst = (double*) erealloc(candGrossEst, numOfCandidates * sizeof(double));
		candGrossSalience = (double*) erealloc(candGrossSalience, numOfCandidates * sizeof(double));
	}
	
	if (estPar->opMode == PRINT_INIT_CANDIDATES) {
		for (i = 0; i < numOfCandidates; i++) {
			printf("%g %g %g\n", time, candGrossEst[i], candGrossSalience[i]);
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* III. Group candidate partials (harmonic series)                                       */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	candHarmFreq = (double**) emalloc(numOfCandidates * sizeof(double*));
	candHarmAmp = (double**) emalloc(numOfCandidates * sizeof(double*));

	maxHarmNumber = round_to_int(peaksFreq[numOfPeaks - 1] / min(candGrossEst, numOfCandidates));

	for (i = 0; i < numOfCandidates; i++) {
		candHarmFreq[i] = (double*) emalloc(maxHarmNumber * sizeof(double));
		candHarmAmp[i] = (double*) emalloc(maxHarmNumber * sizeof(double));
	}

	for (i = 0; i < numOfCandidates; i++) {
		for (j = 0; j < maxHarmNumber; j++) {
			candHarmFreq[i][j] = 0;
			candHarmAmp[i][j] = 0;
		}
	}

	candAverageParcialEnergy = (double*) emalloc(numOfCandidates * sizeof(double));
	candTotalEnergy = (double*) emalloc(numOfCandidates * sizeof(double));
	candFreqError = (double*) emalloc(numOfCandidates * sizeof(double));
	isCandValid = (unsigned char*) emalloc(numOfCandidates * sizeof(unsigned char));
	
	for (i = 0; i < numOfCandidates; i++) {
		for (j = 0; j < numOfPeaks; j++)
		{
			double harm_ratio = peaksFreq[j] / candGrossEst[i];
			unsigned int closest_harm = round_to_int(harm_ratio);
			
			if (absoluteRelativeError(harm_ratio, closest_harm)
				 <=
				estPar->maxRelativeInharmonicity)
			{
				if (peaksAmp[j] > candHarmAmp[i][closest_harm-1]) {
					candHarmFreq[i][closest_harm-1] = peaksFreq[j];
					candHarmAmp[i][closest_harm-1] = peaksAmp[j];
				}
			}
		}
		for (j = 0; j < maxHarmNumber; j++) {
			if (candHarmAmp[i][j] == 0) {
				candHarmFreq[i][j] = candGrossEst[i] * (j+1);
			}
		}

		candTotalEnergy[i] = new_critical_band_sum(candGrossEst[i], candHarmFreq[i], candHarmAmp[i], maxHarmNumber);
	}

	minCandEnergy = estPar->minRelativeEnergy * max(candTotalEnergy, numOfCandidates);

	if (estPar->opMode == PRINT_INIT_CAND_TOTAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			printf("%g %g %g\n", time, candGrossEst[i], candTotalEnergy[i]);
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* IV. Filter out invalid candidates                                                     */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	for (i = 0, numOfValidCandidates = numOfCandidates; i < numOfCandidates; i++) {
		if (candTotalEnergy[i] < minCandEnergy)	{
			isCandValid[i] = NO;

			candTotalEnergy[i] = 0;
			candAverageParcialEnergy[i] = 0;
			
			numOfValidCandidates--;
		}
		else
		{
			isCandValid[i] = YES;
			/*
			candFreqError[i] = f0_error(average_f0_est(candHarmFreq[i], candHarmAmp[i], maxHarmNumber), candHarmFreq[i], candHarmAmp[i], maxHarmNumber);
			*/

			/*
			candAverageParcialEnergy[i] = average2(candHarmAmp[i], maxHarmNumber);
			*/
			/*
			spec_smooth(candHarmFreq[i], candHarmAmp[i], maxHarmNumber);
			*/
			candAverageParcialEnergy[i] = critical_band_average2(candHarmAmp[i], maxHarmNumber);
			candFreqError[i] = totalAbsF0Error(candGrossEst[i], candHarmFreq[i], candHarmAmp[i], candTotalEnergy[i], maxHarmNumber);
		}
	}

	if (numOfValidCandidates == 0) {
		return ERROR;
	}
	
	if (estPar->opMode == PRINT_VALID_CAND_GROSS_FREQ_AND_TOTAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candTotalEnergy[i]);
			}
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	
	if (estPar->opMode == PRINT_VALID_CAND_GROSS_FREQ_AND_AVERAGE_PARTIAL_ENERGY) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candAverageParcialEnergy[i]);
			}
		}
		
		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	if (estPar->opMode == PRINT_VALID_CAND_ERROR) {
		for (i = 0; i < numOfCandidates; i++) {
			if (isCandValid[i] == YES) {
				printf("%g %g %g\n", time, candGrossEst[i], candFreqError[i]);
			}
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* V. Rank candidates                                                                    */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	F0EstIndex = max_index(candAverageParcialEnergy, numOfCandidates);
	
	if (estPar->opMode == PRINT_F0_HARMONIC_SERIES) {
		for (j = 0; j < maxHarmNumber; j++) {
			if (candHarmAmp[F0EstIndex][j]) {
				printf("%g %g %g %u\n", time, candHarmFreq[F0EstIndex][j], candHarmAmp[F0EstIndex][j], j+1);
			}
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
			
	
	if (estPar->opMode == PRINT_GROSS_F0_ESTIMATE)
	{
		printf("%g %g %g\n", time, candGrossEst[F0EstIndex], candAverageParcialEnergy[F0EstIndex]);

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* VI. Refine frequency estimate                                                         */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	/*
	F0EstFreq = average_f0_est(candHarmFreq[F0EstIndex], candHarmAmp[F0EstIndex], maxHarmNumber);
	*/
	F0EstFreq = critical_band_average_f0_estimate(candHarmFreq[F0EstIndex], candHarmAmp[F0EstIndex], maxHarmNumber);


	/* TODO: usar spec_smooth() */

	if (estPar->opMode == PRINT_REFINED_F0_ESTIMATE)
	{
		if (F0EstFreq > 0) {
			printf("%g %g %g\n", time, F0EstFreq, candTotalEnergy[F0EstIndex]);
		} else {
			fprintf(stderr, "%g %g\n", time, F0EstFreq);
		}

		free(candGrossEst);
		free(candGrossSalience);
		free(candAverageParcialEnergy);
		free(candTotalEnergy);
		free(candFreqError);
		free_2d((void**)candHarmFreq, numOfCandidates);
		free_2d((void**)candHarmAmp, numOfCandidates);
		
		return OK;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* VII. Free allocated memory                                                            */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(candGrossEst);
	free(candGrossSalience);
	free(candAverageParcialEnergy);
	free(candTotalEnergy);
	free(candFreqError);
	free_2d((void**)candHarmFreq, numOfCandidates);
	free_2d((void**)candHarmAmp, numOfCandidates);


	return F0EstFreq;
}
