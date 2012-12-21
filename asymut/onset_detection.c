#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include "complex.h"
#include "stdft.h"
#include "config_constants.h"
#include "config_variables.h"
#include "numeric_constants.h"
#include "trigonometric_constants.h"
#include "aux_fun.h"
#include "spec_fun.h"

/********************************************************************************************/
/* LOG_POWER_UNPREDICT                                                                      */
/********************************************************************************************/
double log_power_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
    unsigned int b, num_of_bins = buf->fft_p->num_of_bins;
    unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;
    double total_unpred;

    for (b = 0, total_unpred = 0; b < num_of_bins; b++)
    	total_unpred += fabs(
             				 buf->spec_data[chan][spec_index].log_power.content[b]
             				 -
             				 buf->spec_data[chan][spec_index - 1].log_power.content[b]
                  			);

    return total_unpred;
}

/********************************************************************************************/
/* INCR_LOG_POWER_UNPREDICT                                                                 */
/********************************************************************************************/
/* currently the only one to combine information from critical bands in a special way */
/* TODO: clarify, through the renaming of the variables, the distinction between */
/*       FFT bands (bin) and critical bands */
double incr_log_power_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
    const double critical_bands_per_octave = 3, freq_res = buf->fft_p->freq_resolution;
    double total_increase, curr_mag, prev_mag;
    
    unsigned int k, b, num_of_bins = buf->fft_p->num_of_bins, critical_band,
          			num_of_critical_bands = (critical_bands_per_octave
             								 * log2(buf->stream->sample_rate/freq_res)),
                    band_num_of_bins[num_of_critical_bands];
    unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;

    double curr_band_energy[num_of_critical_bands], prev_band_energy[num_of_critical_bands],
    	   curr_total_energy = buf->spec_data[chan][spec_index].log_power.total_sum,
    	   band_internal_movement[num_of_critical_bands];

	/*
    if (curr_total_energy <= 0)
    	return 0;
	*/

    for (b = 0; b < num_of_critical_bands; b++)
    {
    	curr_band_energy[b] = 0;
    	prev_band_energy[b] = 0;
    	band_num_of_bins[b] = 0;
    	band_internal_movement[b] = 0;
    }
    
    const unsigned int first_bin = ceil(MIN_AUDIBLE_FREQUENCY/freq_res) - 1,
    					last_bin = ceil(MAX_AUDIBLE_FREQUENCY/freq_res) - 1;

    for (k = first_bin; k <= last_bin; k++)
    {
        curr_mag = buf->spec_data[chan][spec_index].power.content[k]; /* TEST */
        prev_mag = buf->spec_data[chan][spec_index - 1].power.content[k]; /* TEST */
        critical_band = critical_bands_per_octave * log2(freq_res*(k+1)/MIN_AUDIBLE_FREQUENCY);

       	band_num_of_bins[critical_band]++;

       	band_internal_movement[critical_band] += fabs(curr_mag - prev_mag);
       	
        curr_band_energy[critical_band] += curr_mag;
        prev_band_energy[critical_band] += prev_mag;
   	}
   	
   	/*
   	double curr_time = ((double)buf->first_sample_file_pos + (unpred_index+2)*buf->fft_p->step) /buf->stream->sample_rate;
   	*/
   	for (b = 0; b < num_of_critical_bands; b++)
        if (band_num_of_bins[b] >= 1)
        {
        	/*
        	band_internal_movement[b] /= fabs(curr_band_energy[b] - prev_band_energy[b]);
        	*/

    	    curr_band_energy[b] = log10((double)curr_band_energy[b]/(double)band_num_of_bins[b]);
    		prev_band_energy[b] = log10((double)prev_band_energy[b]/(double)band_num_of_bins[b]);
    		
    		/*
    		printf("%g %u %g\n", curr_time, b, curr_band_energy[b] - prev_band_energy[b]);
    		*/
      	}
	/*
	printf("\n\n");
	return 0;
	*/

    for (b = 0, total_increase = 0; b < num_of_critical_bands; b++) {
    	if (curr_band_energy[b] > prev_band_energy[b])
    	{
         /*
        	total_increase += sqrt( pow( (curr_band_energy[b] - prev_band_energy[b]), 2) +
         						pow(band_internal_movement[b], 2) );
         						*/
			/*
         	total_increase += (curr_band_energy[b] - prev_band_energy[b]);
         	*/
   		}
   		if (band_num_of_bins[b] >= 1) {
	       	total_increase += curr_band_energy[b];
   		}
    }

    return total_increase; /* * curr_total_energy;*/
}

/********************************************************************************************/
/* MKL_UNPREDICT                                                                            */
/********************************************************************************************/
/* based on (Stephen Hainsworth, 2003) paper */
/* suggested parameters: 4096 samples window, 512 samples step */
double mkl_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned int spec_index = buf->unpredictability_to_spec_offset + unpred_index;
	const double lowest_freq = 30, highest_freq = 5e3;
	
	unsigned int k, lowest_bin, highest_bin;
	double total_distance;
	
	lowest_bin = ceil(lowest_freq/buf->fft_p->freq_resolution)-1;
	highest_bin = floor(highest_freq/buf->fft_p->freq_resolution)-1;
	
	for(k = lowest_bin, total_distance = 0; k <= highest_bin; k++) {
		double bin_distance = log10(buf->spec_data[chan][spec_index].magnitude.content[k]
								/buf->spec_data[chan][spec_index-1].magnitude.content[k])/LOGARITHM_OF_2_IN_BASE_10;
		
		if (bin_distance > 0) {
			total_distance += bin_distance;
		}
	}
	return total_distance;
}

/********************************************************************************************/
/* SPECTRAL_DIFFERENCE_UNPREDICT                                                            */
/********************************************************************************************/
/* based on (Juan Bello, 2005) paper */
/* suggested parameters: 4096 samples window, 512 samples step */
double spectral_difference_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned int spec_index = buf->unpredictability_to_spec_offset + unpred_index,
						num_of_bins = buf->fft_p->num_of_bins;
	
	unsigned int k;
	double total_distance;
	
	for(k = total_distance = 0; k < num_of_bins; k++) {
		double bin_distance = buf->spec_data[chan][spec_index].magnitude.content[k]
								- buf->spec_data[chan][spec_index-1].magnitude.content[k];
		
		if (bin_distance > 0) {
			total_distance += (bin_distance*bin_distance);
		}
	}
	return total_distance;
}

/********************************************************************************************/
/* ERB_FRAMEWORK_UNPREDICT                                                                  */
/********************************************************************************************/
/* based on (Nick Collins, 2005) paper
 * ERB stands for Equivalent Rectangular Bandwidth scale*/
double erb_framework_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index,
		unsigned char mode)
{
	const unsigned char num_of_bands = 40; /* number of ERB bands */
    const unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;
    const double freq_res = buf->fft_p->freq_resolution;

    unsigned int k, b, band_num_of_bins[num_of_bands];

    double curr_band_energy[num_of_bands], prev1_band_energy[num_of_bands],
    		prev2_band_energy[num_of_bands], total_increase;

    for (b = 0; b < num_of_bands; b++)
    {
    	band_num_of_bins[b] = 0;
    	curr_band_energy[b] = 0;
    	prev1_band_energy[b] = 0;
    	prev2_band_energy[b] = 0;
    }
    
    const unsigned int first_bin = ceil(MIN_AUDIBLE_FREQUENCY/freq_res) - 1,
    					last_bin = floor(MAX_AUDIBLE_FREQUENCY/freq_res) - 1;

    for (k = first_bin; k <= last_bin; k++)
    {
        b = floor(21.4 * log10(4.37e-3*(k+1)*freq_res + 1) + 0.5);
        
        if (b >= num_of_bands) break;

       	band_num_of_bins[b]++;
       	
        curr_band_energy[b] += buf->spec_data[chan][spec_index].power.content[k];
        prev1_band_energy[b] += buf->spec_data[chan][spec_index-1].power.content[k];
        prev2_band_energy[b] += buf->spec_data[chan][spec_index-2].power.content[k];
   	}
   	
   	for (b = total_increase = 0; b < num_of_bands; b++)
        if (band_num_of_bins[b] >= 1)
        {
    	    if (mode == 0) {
	    	    curr_band_energy[b] = log10((double)curr_band_energy[b]/band_num_of_bins[b]);
	    	    
	    		prev1_band_energy[b] = log10((double)prev1_band_energy[b]/band_num_of_bins[b]);
	    		prev2_band_energy[b] = log10((double)prev2_band_energy[b]/band_num_of_bins[b]);

		    	double band_increase = curr_band_energy[b] - 0.5*(prev1_band_energy[b] + prev2_band_energy[b]);
		    	
		    	if (band_increase > 0) total_increase += band_increase;
    	    } else {
    	    	total_increase += log10((double)curr_band_energy[b]/band_num_of_bins[b]);
    	    }
      	}

    return 10*total_increase; /* Bel to decibel */
}

/********************************************************************************************/
/* ERB_UNPREDICT                                                                            */
/********************************************************************************************/
/* based on (Nick Collins, 2005) paper
 * ERB stands for Equivalent Rectangular Bandwidth scale*/
double old_erb_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	return erb_framework_unpredict(buf, chan, unpred_index, 1);	
}

/********************************************************************************************/
/* BARK_UNPREDICT                                                                           */
/********************************************************************************************/
/* based on (Nick Collins, 2005) paper,
 * with Bark instead of Equivalent Rectangular Bandwidth scale */
double bark_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned char num_of_bands = 26; /* number of Bark bands */
    const unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;
    const double freq_res = buf->fft_p->freq_resolution;

    unsigned int k, b, band_num_of_bins[num_of_bands];

    double curr_band_energy[num_of_bands], prev1_band_energy[num_of_bands],
    		prev2_band_energy[num_of_bands], total_increase;

    for (b = 0; b < num_of_bands; b++)
    {
    	band_num_of_bins[b] = 0;
    	curr_band_energy[b] = 0;
    	prev1_band_energy[b] = 0;
    	prev2_band_energy[b] = 0;
    }
    
    const unsigned int first_bin = ceil(MIN_AUDIBLE_FREQUENCY/freq_res) - 1,
    					last_bin = floor(MAX_AUDIBLE_FREQUENCY/freq_res) - 1;

    for (k = first_bin; k <= last_bin; k++)
    {
        b = floor(hertz_to_bark((k+1)*freq_res) + 0.5);
        
        if (b >= num_of_bands) break;

       	band_num_of_bins[b]++;
       	
        curr_band_energy[b] += buf->spec_data[chan][spec_index].power.content[k];
        prev1_band_energy[b] += buf->spec_data[chan][spec_index-1].power.content[k];
        prev2_band_energy[b] += buf->spec_data[chan][spec_index-2].power.content[k];
   	}
   	
   	for (b = total_increase = 0; b < num_of_bands; b++)
        if (band_num_of_bins[b] >= 1)
        {
    	    curr_band_energy[b] = log10((double)curr_band_energy[b]/band_num_of_bins[b]);
    		prev1_band_energy[b] = log10((double)prev1_band_energy[b]/band_num_of_bins[b]);
    		prev2_band_energy[b] = log10((double)prev2_band_energy[b]/band_num_of_bins[b]);

	    	double band_increase = curr_band_energy[b] - 0.5*(prev1_band_energy[b] + prev2_band_energy[b]);
	    	
	    	if (band_increase > 0) total_increase += 10*band_increase; /* Bel to decibel */
      	}

    return total_increase;
}

/********************************************************************************************/
/* RMS_UNPREDICT                                                                            */
/********************************************************************************************/
double rms_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned int num_of_bins = buf->fft_p->num_of_bins,
							spec_index = unpred_index + buf->unpredictability_to_spec_offset;
	unsigned int i;
	double curr_energy, prev_energy;
	
	for (curr_energy = prev_energy = 0, i = 0; i < num_of_bins; i++) {
		prev_energy += buf->spec_data[chan][spec_index-1].power.content[i];
		curr_energy += buf->spec_data[chan][spec_index].power.content[i];
	}
	if (curr_energy > prev_energy) {
		return sqrt(curr_energy) - sqrt(prev_energy);
	}
	return 0;
}

/********************************************************************************************/
/* WEIGHTED_RMS_UNPREDICT                                                                   */
/********************************************************************************************/
double weighted_rms_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned int num_of_bins = buf->fft_p->num_of_bins;
	unsigned int k;
	double energy1, energy2;
	
	for (energy1 = energy2 = 0, k = 0; k <= num_of_bins; k++) {
		energy1 += buf->spec_data[chan][unpred_index+1].power.content[k] * (k+1);
		energy2 += buf->spec_data[chan][unpred_index+2].power.content[k] * (k+1);
	}
	
	if (energy2 > energy1) return sqrt(energy2)-sqrt(energy1);
	return 0;
}

/********************************************************************************************/
/* NEW_UNPREDICT                                                                            */
/********************************************************************************************/
double new_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
    const unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;
    const double critical_bands_per_octave = 3, freq_res = buf->fft_p->freq_resolution;

    double band_lp_dif, total_increase, mag2, mag3, lp2, lp3, uw1, uw2, uw3, cum1, cum2, cum3;

    unsigned int b, num_of_bins = buf->fft_p->num_of_bins, critical_band,
          			num_of_critical_bands = critical_bands_per_octave
             								 * log2(num_of_bins * freq_res / 20),
                    band_num_of_bins[num_of_critical_bands];

    double curr_band_mag_sum[num_of_critical_bands], prev_band_mag_sum[num_of_critical_bands],
           curr_band_lp_sum[num_of_critical_bands], prev_band_lp_sum[num_of_critical_bands],
    	   curr_total_energy = buf->spec_data[chan][spec_index].log_power.total_sum,
    	   band_internal_movement[num_of_critical_bands],
           band_ph_unpred[num_of_critical_bands],
           x[num_of_critical_bands],
           y[num_of_critical_bands],
           z[num_of_critical_bands],
           sam[num_of_critical_bands];

	signed long n1 = buf->first_sample_file_pos + buf->unpredictability_to_sample_offset
					+ buf->fft_p->step * (spec_index - 2),
				n2 = buf->first_sample_file_pos + buf->unpredictability_to_sample_offset
 					+ buf->fft_p->step * (spec_index - 1),
 				n3 = buf->first_sample_file_pos + buf->unpredictability_to_sample_offset
 					+ buf->fft_p->step * spec_index;

    if (curr_total_energy <= 0)
    	return 0;

    for (b = 0; b < num_of_critical_bands; b++)
    {
     	prev_band_mag_sum[b] = 0;
    	prev_band_lp_sum[b] = 0;

    	curr_band_mag_sum[b] = 0;
    	curr_band_lp_sum[b] = 0;

    	band_num_of_bins[b] = 0;
    	band_internal_movement[b] = 0;
    	band_ph_unpred[b] = 0;
    	
    	x[b] = 0;
    	y[b] = 0;
    	z[b] = 0;
    	sam[b] = 0;
    }

    for (b = 20/freq_res; b < 20000/freq_res; b++)
    {
        cum1 = (TWO_PI * (b + 1) * n1) / buf->fft_p->window.length;
     	uw1 = buf->spec_data[chan][spec_index - 2].phase.content[b] + cum1;

     	cum2 = (TWO_PI * (b + 1) * n2) / buf->fft_p->window.length;
        uw2 = buf->spec_data[chan][spec_index - 1].phase.content[b] + cum2;
        lp2 = buf->spec_data[chan][spec_index - 1].log_power.content[b];
        mag2 = pow(10, lp2);

     	cum3 = (TWO_PI * (b + 1) * n3) / buf->fft_p->window.length;
        uw3 = buf->spec_data[chan][spec_index].phase.content[b] + cum3;
        lp3 = buf->spec_data[chan][spec_index].log_power.content[b];
        mag3 = pow(10, lp3);

        critical_band = critical_bands_per_octave * log2(freq_res * (b+1) / 20);

       	band_num_of_bins[critical_band]++;

       	band_ph_unpred[critical_band] += fabs(princ_det(uw3 - 2*uw2 + uw1))* lp3;

       	band_internal_movement[critical_band] += fabs(lp3 - lp2);
        prev_band_lp_sum[critical_band] += lp2;
       	curr_band_lp_sum[critical_band] += lp3;

        curr_band_mag_sum[critical_band] += mag3;
        prev_band_mag_sum[critical_band] += mag2;

        x[critical_band] += (mag2 * mag3);
        y[critical_band] += (mag2 * mag2);
        z[critical_band] += (mag3 * mag3);
   	}

   	for (b = 0; b < num_of_critical_bands; b++)
        if (band_num_of_bins[b] >= 1)
        {
        	band_ph_unpred[b] /= curr_band_lp_sum[b];

        	if (curr_band_mag_sum[b] > 1)
        		curr_band_lp_sum[b] = log10(curr_band_mag_sum[b]);
        	else
        		curr_band_lp_sum[b] = 0;

        	if (prev_band_mag_sum[b] > 1)
           		prev_band_lp_sum[b] = log10(prev_band_mag_sum[b]);
           	else
           		prev_band_lp_sum[b] = 0;

    		/*
    		band_internal_movement[b] /= band_num_of_bins[b];
    		*/
    		
    		if (y[b] > 0 && z[b] > 0)
    		{
    			sam[b] = acos((x[b] / (sqrt(y[b]) * sqrt(z[b])))
    					/ ((y[b] > z[b] ? y[b] : z[b]) / (sqrt(y[b]) * sqrt(z[b]))));
    	    }
    		else
    		{
    			if (y[b] > 0 || z[b] > 0)
    				sam[b] = PI/2;
    			else
    				sam[b] = 0;
      		}
      	}

    for (b = 0, total_increase = 0; b < num_of_critical_bands; b++)
    {
    	if (curr_band_mag_sum[b] > prev_band_mag_sum[b])
    	{
         	band_lp_dif = curr_band_lp_sum[b] - prev_band_lp_sum[b];
         	/*
          	if (band_lp_dif < 1)
         		band_lp_dif = 1;

    		if (band_internal_movement[b] < 1)
    			band_internal_movement[b] = 1;

    		if (curr_band_lp_sum[b] < 1)
    			curr_band_lp_sum[b] = 1;

       		if (band_ph_unpred[b] < 1)
       			band_ph_unpred[b] = 1;
       			*/
       	       			/*
            total_increase += band_lp_dif
            				   * band_internal_movement[b]
                   				* curr_band_lp_sum[b]
            				   * band_ph_unpred[b];
            				   */
            total_increase += band_lp_dif
                   			  * curr_band_lp_sum[b]
            				  * band_ph_unpred[b]
                  			  * pow(sam[b], 16)
                          ;
        }
        else
        {
            /*
        	aux = (curr_band_lp_sum[b] / prev_band_lp_sum[b]) * band_internal_movement[b]
                   * curr_band_lp_sum[b] * band_ph_unpred[b];

            if (aux > 0)
            	total_increase += aux;
            band_lp_dif = pow(pow(sam[b], 16) * curr_band_lp_sum[b], 2);
            if (band_lp_dif != 0)
            	printf("Ah! %g\n", band_lp_dif);
            total_increase += band_lp_dif;
            */
        }
    }

    return total_increase;
}


/********************************************************************************************/
/* PHASE_UNPREDICT                                                                          */
/********************************************************************************************/
double phase_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	unsigned int b;
  	double uw, uw1, uw2, total_deviation;
  	
	for (b = 0, total_deviation = 0; b < buf->fft_p->num_of_bins; b++)
	{
		double phase_dif;

     	uw = buf->spec_data[ch][index].unwrapped_phase.content[b];
        uw1 = buf->spec_data[ch][index + 1].unwrapped_phase.content[b];
        uw2 = buf->spec_data[ch][index + 2].unwrapped_phase.content[b];

		phase_dif = fmod((2*uw1 - uw), TWO_PI) - fmod(uw2, TWO_PI);
		
       	total_deviation += fabs(phase_dif);
	}

	return total_deviation / buf->fft_p->num_of_bins;
}

/********************************************************************************************/
/* MODIFIED_PHASE_UNPREDICT                                                                 */
/********************************************************************************************/
double modified_phase_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	unsigned int b;

  	double uw, uw1, uw2, p1, p2, total_deviation;
  	
	for (b = 0, total_deviation = 0; b < buf->fft_p->num_of_bins; b++)
	{
        p1 = buf->spec_data[ch][index + 1].power.content[b];
        p2 = buf->spec_data[ch][index + 2].power.content[b];
        
		if (p2>p1) {
			double phase_dif;

	     	uw = buf->spec_data[ch][index].unwrapped_phase.content[b];
	        uw1 = buf->spec_data[ch][index + 1].unwrapped_phase.content[b];
	        uw2 = buf->spec_data[ch][index + 2].unwrapped_phase.content[b];

			phase_dif = fmod((2*uw1 - uw), TWO_PI) - fmod(uw2, TWO_PI);
			
	       	total_deviation += fabs(phase_dif)*p2;
		}
	}

	return total_deviation / buf->fft_p->num_of_bins;
}

/********************************************************************************************/
/* COMPLEX_UNPREDICT                                                                        */
/********************************************************************************************/
double complex_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	unsigned int b, num_of_bins = buf->fft_p->num_of_bins;
	
  	double total_deviation, target_mag, measured_mag, phase_dif;

    double *spec_uph_content = buf->spec_data[ch][index].unwrapped_phase.content,
    		*spec1_uph_content = buf->spec_data[ch][index + 1].unwrapped_phase.content,
        	*spec2_uph_content = buf->spec_data[ch][index + 2].unwrapped_phase.content,
    		*spec1_mag_content = buf->spec_data[ch][index + 1].magnitude.content,
    		*spec2_mag_content = buf->spec_data[ch][index + 2].magnitude.content;

	for (b = total_deviation = 0; b < num_of_bins; b++)
	{
		target_mag = spec1_mag_content[b];
		measured_mag = spec2_mag_content[b];
		
		phase_dif = fmod((2*spec1_uph_content[b] - spec_uph_content[b]), TWO_PI) - fmod(spec2_uph_content[b], TWO_PI);
		
		total_deviation += sqrt(target_mag*target_mag + measured_mag*measured_mag
								-2*target_mag*measured_mag*cos(phase_dif));
	}

	return total_deviation;
}


/********************************************************************************************/
/* NON_DECREASING_COMPLEX_UNPREDICT                                                         */
/********************************************************************************************/
double non_decreasing_complex_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	unsigned int b, num_of_bins = buf->fft_p->num_of_bins;
	
  	double total_deviation, target_mag, measured_mag, phase_dif;

    double *spec_uph_content = buf->spec_data[ch][index].unwrapped_phase.content,
    		*spec1_uph_content = buf->spec_data[ch][index + 1].unwrapped_phase.content,
        	*spec2_uph_content = buf->spec_data[ch][index + 2].unwrapped_phase.content,
    		*spec1_mag_content = buf->spec_data[ch][index + 1].magnitude.content,
    		*spec2_mag_content = buf->spec_data[ch][index + 2].magnitude.content;

	for (b = total_deviation = 0; b < num_of_bins; b++)
	{
		target_mag = spec1_mag_content[b];
		measured_mag = spec2_mag_content[b];
		
		if (measured_mag >= target_mag) {
			phase_dif = fmod((2*spec1_uph_content[b] - spec_uph_content[b]), TWO_PI) - fmod(spec2_uph_content[b], TWO_PI);
			
			total_deviation += sqrt(target_mag*target_mag + measured_mag*measured_mag
									-2*target_mag*measured_mag*cos(phase_dif));
		}
	}

	return total_deviation;
}


/********************************************************************************************/
/* OLD_SLOW_COMPLEX_UNPREDICT                                                               */
/********************************************************************************************/
double old_slow_complex_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	unsigned int b, num_of_bins = buf->fft_p->num_of_bins;
	
  	double target_ph, measured_ph, target_mag, measured_mag, total_deviation, real_dif, im_dif;

    double *spec_uph_content = buf->spec_data[ch][index].unwrapped_phase.content,
    		*spec1_uph_content = buf->spec_data[ch][index + 1].unwrapped_phase.content,
        	*spec2_uph_content = buf->spec_data[ch][index + 2].unwrapped_phase.content,
    		*spec1_mag_content = buf->spec_data[ch][index + 1].magnitude.content,
    		*spec2_mag_content = buf->spec_data[ch][index + 2].magnitude.content;

	for (b = total_deviation = 0; b < num_of_bins; b++)
	{
     	target_mag = spec1_mag_content[b];
        measured_mag = spec2_mag_content[b];
        
        target_ph = fmod((2*spec1_uph_content[b] - spec_uph_content[b]), TWO_PI);
        measured_ph = fmod(spec2_uph_content[b], TWO_PI);
        
        real_dif = target_mag * cos(target_ph) - measured_mag * cos(measured_ph);
        im_dif = target_mag * sin(target_ph) - measured_mag * sin(measured_ph);

        total_deviation += sqrt(real_dif*real_dif + im_dif*im_dif);
	}

	return total_deviation;
}

/********************************************************************************************/
/* PSYCHOACOUSTIC_UNPREDICT                                                                 */
/********************************************************************************************/
double psychoacoustic_unpredict(Buffer* buf, unsigned char ch, unsigned int index)
{
	/* TODO: implement something more realistic than ILP */
	
	
	return 0;
}



/********************************************************************************************/
/* ERB_UNPREDICT                                                                  */
/********************************************************************************************/
/* based on (Nick Collins, 2005) paper
 * ERB stands for Equivalent Rectangular Bandwidth scale*/
double erb_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index)
{
	const unsigned char num_of_bands = 40; /* number of ERB bands */
    const unsigned int spec_index = unpred_index + buf->unpredictability_to_spec_offset;
    const double freq_res = buf->fft_p->freq_resolution;

    unsigned int k, b, band_num_of_bins[num_of_bands];

    double curr_band_energy[num_of_bands], prev1_band_energy[num_of_bands],
    		prev2_band_energy[num_of_bands], total_increase;

    for (b = 0; b < num_of_bands; b++)
    {
    	band_num_of_bins[b] = 0;
    	curr_band_energy[b] = 0;
    	prev1_band_energy[b] = 0;
    	prev2_band_energy[b] = 0;
    }
    
    const unsigned int first_bin = ceil(MIN_AUDIBLE_FREQUENCY/freq_res) - 1,
    					last_bin = floor(MAX_AUDIBLE_FREQUENCY/freq_res) - 1;

    for (k = first_bin; k <= last_bin; k++)
    {
        b = floor(21.4 * log10(4.37e-3*(k+1)*freq_res + 1) + 0.5);
        
        if (b >= num_of_bands) break;

       	band_num_of_bins[b]++;
       	
        curr_band_energy[b] += buf->spec_data[chan][spec_index].power.content[k];
        prev1_band_energy[b] += buf->spec_data[chan][spec_index-1].power.content[k];
        prev2_band_energy[b] += buf->spec_data[chan][spec_index-2].power.content[k];
   	}
   	
   	for (b = total_increase = 0; b < num_of_bands; b++)
        if (band_num_of_bins[b] >= 1)
        {
    	    if (0 == 0) {
	    	    curr_band_energy[b] = log10((double)curr_band_energy[b]/band_num_of_bins[b]);
	    	    
	    		prev1_band_energy[b] = log10((double)prev1_band_energy[b]/band_num_of_bins[b]);
	    		prev2_band_energy[b] = log10((double)prev2_band_energy[b]/band_num_of_bins[b]);

		    	double band_increase = curr_band_energy[b] - 0.5*(prev1_band_energy[b] + prev2_band_energy[b]);
		    	
		    	if (band_increase > 0) total_increase += band_increase;
    	    } else {
    	    	total_increase += log10((double)curr_band_energy[b]/band_num_of_bins[b]);
    	    }
      	}

    return 10*total_increase; /* Bel to decibel */
}
