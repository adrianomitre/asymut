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

static double MAG_NORM_FACTOR, POW_NORM_FACTOR, LP_NORM_PARCEL;

/********************************************************************************************/
/* INIT_NORMALIZATION                                                                       */
/********************************************************************************************/
void init_normalization(const Buffer *buf)
{
	if (NORMALIZATION_ENABLED == FALSE) {
		LP_NORM_PARCEL = 0;
		POW_NORM_FACTOR = 1;
		MAG_NORM_FACTOR = 1;

		return;
	}

	const double max_amplitude = buf->fft_p->window.total_area * 0.5
									* (pow(2, buf->stream->bit_depth - 1) + 0.5),
					max_power = max_amplitude * max_amplitude,
   					max_log_power = 10 * log10(max_power);

	LP_NORM_PARCEL = NORMALIZATION_LEVEL - max_log_power;
	POW_NORM_FACTOR = 1/max_power;
	MAG_NORM_FACTOR = 1/max_amplitude;
}


/********************************************************************************************/
/* COPY_SPECTRA                                                                             */
/********************************************************************************************/
void copy_spectra(Buffer* buf, const unsigned char chan, const unsigned int src,
																	const unsigned int dest)
{
    const unsigned int num_of_bins = buf->fft_p->num_of_bins;

    copy_spectrum(&buf->spec_data[chan][src].phase,
    								&buf->spec_data[chan][dest].phase, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].unwrapped_phase,
    								&buf->spec_data[chan][dest].unwrapped_phase, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].magnitude,
    								&buf->spec_data[chan][dest].magnitude, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].power,
    								&buf->spec_data[chan][dest].power, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].log_power,
    								&buf->spec_data[chan][dest].log_power, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].warped_denoised,
									&buf->spec_data[chan][dest].warped_denoised, num_of_bins);
    copy_spectrum(&buf->spec_data[chan][src].power_noise,
    								&buf->spec_data[chan][dest].power_noise, num_of_bins);
}

/********************************************************************************************/
/* COPY_SPECTRUM                                                                            */
/********************************************************************************************/
void copy_spectrum(Spectrum *src, Spectrum *dest, const unsigned int num_of_bins)
{
    memmove(dest->content, src->content, num_of_bins * sizeof(double));
    dest->max_index = src->max_index;
    dest->total_sum = src->total_sum;
}

/********************************************************************************************/
/* MAG_WARP                                                                                 */
/********************************************************************************************/
double mag_warp(double original[], double warped[], const unsigned int min,
											const unsigned int max, const unsigned int n)
{
    unsigned int k;
    double sum;

    if (min > max || max >= n)
    {
        fprintf(stderr, "Fatal error in mag_warp: inconsistent min, max, "
        				"n values (%u, %u, %u). Exiting...\n", min, max, n);
        exit(-1);
   	}

	for (sum = 0.0, k = min; k <= max; k++)
    	sum += pow(original[k], 0.33333333333333333333333333333333);
    sum = pow(sum/(max - min + 1), 3);

    for (k = 0; k < n; k++)
    {
        if (sum > 0) warped[k] = log(1 + original[k]/sum);
        else warped[k] = 0;
    }

    return sum;
}

/********************************************************************************************/
/* NOISE_ESTIMATE_AND_SUPRESS                                                               */
/********************************************************************************************/
void noise_estimate_and_supress(double spec[], double power_noise[], double warp_factor,
																	Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands, last_band_index = num_of_bands - 1;
    unsigned int k, num_of_bins = aud_mod->fft_p->num_of_bins, first_interp, last_interp,
    			   last_bin;
    double band_mean[num_of_bands], aux_sum, freq_res = aud_mod->fft_p->freq_resolution;

    for (b = 0; b < num_of_bands; b++)
   	{
   		if (aud_mod->band[b].num_of_bins != 0)
   		{
	        for (k = aud_mod->band[b].first_bin, aux_sum = 0.0,
	        	 last_bin = aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins - 1;
	        	 k <= last_bin;
	        	 k++)
	      	 {
	        	 aux_sum += (spec[k] * aud_mod->band[b].response[k]);
	      	 }
	      	 band_mean[b] = aux_sum / aud_mod->band[b].response_total_sum;
   		}
		else band_mean[b] = 0;
    }

    first_interp = ceil(aud_mod->band[0].center / freq_res) - 1;
    last_interp = floor(aud_mod->band[last_band_index].center / freq_res) - 1;

    /* in reality, the next 'for' may be omitted for the sake of speed
    	as the extremes of the spectrum will not be used for pitch-estimating */
    for (k = 0; k < first_interp; k++)
    {
        power_noise[k] = band_mean[0];

    	spec[k] -= power_noise[k];

     	if (spec[k] < 0)
        	spec[k] = 0;

        power_noise[k] = (exp(power_noise[k]) - 1) * warp_factor;
    }

    for (k = first_interp, b = 0; k <= last_interp; k++)
    {
        while (b < (num_of_bands - 1) && (k + 1) * freq_res > aud_mod->band[b + 1].center)
        {
        	b++;
       	}

        power_noise[k] = ( band_mean[b] + (band_mean[b + 1] - band_mean[b])
        			* ( ((k + 1) * freq_res - aud_mod->band[b].center)
           			/ (aud_mod->band[b + 1].center - aud_mod->band[b].center) ) );

        spec[k] -= power_noise[k];

        if (spec[k] < 0)
        	spec[k] = 0;

        power_noise[k] = (exp(power_noise[k]) - 1) * warp_factor;
    }

    /* in reality, the next 'for' may be omitted for the sake of speed
    	as the extremes of the spectrum will not be used for pitch-estimating */
    for (k = last_interp + 1; k < num_of_bins; k++)
    {
    	power_noise[k] = band_mean[last_band_index];

    	spec[k] -= power_noise[k];

    	if (spec[k] < 0)
        	spec[k] = 0;

        power_noise[k] = (exp(power_noise[k]) - 1) * warp_factor;
    }
}

/********************************************************************************************/
/* MAG_WARP_DENOISE                                                                         */
/********************************************************************************************/
void mag_warp_denoise(double power_spec[], double warped_denoised[], double power_noise[],
																	Auditory_Model* aud_mod)
{
    const unsigned int min = aud_mod->band[0].first_bin,
    					 max = aud_mod->band[aud_mod->num_of_bands - 1].first_bin
          							+ aud_mod->band[aud_mod->num_of_bands - 1].num_of_bins;
    double warp_factor;

	warp_factor = mag_warp(power_spec, warped_denoised, min, max, aud_mod->fft_p->num_of_bins);
 	noise_estimate_and_supress(warped_denoised, power_noise, warp_factor, aud_mod);

	return;
}

/********************************************************************************************/
/* FUTURELY_FASTER_FFT_PH_UPH_MAG_POW_LP_WD_PN                                                              */
/********************************************************************************************/
void futurely_faster_fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod)
{
	static unsigned char first_call = TRUE;
	static Complex *c_in = NULL, *c_out = NULL;

	unsigned int i, i_plus_one, core_length, total_length, last_bin;
	float *signal;
	double aux, cum, cum_incr, *fft_p_window_content,
			*spec_phase_content, *spec_unwrapped_phase_content, *spec_magnitude_content,
   			*spec_power_content, *spec_log_power_content,
   			out_re, out_im;

	FFT_properties *fft_p;
	Spectra *spec;

	if (buf == NULL || aud_mod == NULL)
	{
        if (!first_call)
        {
         free(c_in);
         free(c_out);
        }
        return;
    }

    signal = &(buf->sample_data[channel][sample_index]);
	fft_p = buf->fft_p;
	fft_p_window_content = &(fft_p->window.content[0]);
	spec = &(buf->spec_data[channel][spec_index]);
	core_length = fft_p->window.length / 2;
	total_length = core_length * fft_p->zero_padding_ratio / 2;
	last_bin = fft_p->num_of_bins - 1;

	if (first_call)
 	{
      c_in = (Complex*) emalloc(total_length * sizeof(Complex));
	  c_out = (Complex*) emalloc(total_length * sizeof(Complex));

	  for (i = core_length; i < total_length; i++)
	  	c_in[i].re = c_in[i].im = 0;

	  first_call = FALSE;
	}

    if (buf->fft_p->window.apodization_function.type == RECTANGULAR)
    {
    	for (i = 0; i < core_length; i++)
    	{
    		c_in[i].re = signal[2*i];
    		c_in[i].im = signal[2*i+1];
    	}
    }
    else
    {
    	for (i = 0; i < core_length; i++)
    	{
    		c_in[i].re = signal[2*i] * fft_p_window_content[2*i];
    		c_in[i].im = signal[2*i+1] * fft_p_window_content[2*i+1];
    	}
    }

    fft(c_in, c_out, total_length);

    spec_phase_content = spec->phase.content;
    spec_magnitude_content = spec->magnitude.content;
    spec_power_content = spec->power.content;
   	spec_log_power_content = spec->log_power.content;
   	spec_unwrapped_phase_content = spec->unwrapped_phase.content;

  	for (i = 0, i_plus_one = 1, cum = cum_incr = (buf->first_sample_file_pos + sample_index)
  																	* TWO_PI / core_length;
   		 i < last_bin;
       	 i++, i_plus_one++, cum += cum_incr)
    {
    	Complex a, b, c, half_i;
    	half_i.re = 0;
    	half_i.im = 0.5;

    	c_copy(c_mult_by_real(c_add(c_out[i_plus_one], c_conjugate(c_out[total_length - i_plus_one])), 0.5), &a);

    	c_copy(c_mult(c_subtract(c_out[i_plus_one], c_conjugate(c_out[total_length - i_plus_one])), half_i), &b);
    	c_copy(c_mult(c_exp(TWO_PI*i_plus_one/total_length), b), &b);

    	c_copy(c_subtract(a, b), &c);

    	out_re =  c.re;
    	out_im = c.im;

        spec_phase_content[i] = atan2(out_im, out_re);
        spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

		aux = out_re * out_re + out_im * out_im;

   		spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   		spec_power_content[i] = aux * POW_NORM_FACTOR;
   		spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;
    }
    /* deals with the problematic last bin separately */
	out_re =  c_out[i_plus_one].re;
	out_im = c_out[i_plus_one].im;

	spec_phase_content[i] = atan2(out_im, out_re);
    spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

	aux = (out_re * out_re + out_im * out_im)/4;


	spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   	spec_power_content[i] = aux * POW_NORM_FACTOR;
   	spec_log_power_content[i] = log10(aux) + LP_NORM_PARCEL;

   	mag_warp_denoise(spec->power.content, spec->warped_denoised.content,
    													spec->power_noise.content, aud_mod);
    calc_spectra_extra_data(spec, last_bin + 1);
}


/********************************************************************************************/
/* FFT_PH_UPH_MAG_POW_LP_WD_PN                                                              */
/********************************************************************************************/
void fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod)
{
	static unsigned char first_call = TRUE;
	static Complex *c_in = NULL, *c_out = NULL;

	unsigned int i, i_plus_one, core_length, total_length, last_bin;
	float *signal;
	double aux, cum, cum_incr, *fft_p_window_content,
			*spec_phase_content, *spec_unwrapped_phase_content, *spec_magnitude_content,
   			*spec_power_content, *spec_log_power_content;

	FFT_properties *fft_p;
	Spectra *spec;

	if (buf == NULL || aud_mod == NULL)
	{
        if (!first_call)
        {
         free(c_in);
         free(c_out);
        }
        return;
    }

    signal = &(buf->sample_data[channel][sample_index]);
	fft_p = buf->fft_p;
	fft_p_window_content = &(fft_p->window.content[0]);
	spec = &(buf->spec_data[channel][spec_index]);
	core_length = fft_p->window.length;
	total_length = core_length * fft_p->zero_padding_ratio;
	last_bin = fft_p->num_of_bins - 1;

	if (first_call)
 	{
      c_in = (Complex*) emalloc(total_length * sizeof(Complex));
	  c_out = (Complex*) emalloc(total_length * sizeof(Complex));

   	  for (i = 0; i < core_length; i++)
    	c_in[i].im = 0;

	  for (i = core_length; i < total_length; i++)
	  	c_in[i].re = c_in[i].im = 0;

	  first_call = FALSE;
	}

   	for (i = 0; i < core_length; i++)
    		c_in[i].re = signal[i] * fft_p_window_content[i];

    fft(c_in, c_out, total_length);

    spec_phase_content = spec->phase.content;
    spec_magnitude_content = spec->magnitude.content;
    spec_power_content = spec->power.content;
   	spec_log_power_content = spec->log_power.content;
   	spec_unwrapped_phase_content = spec->unwrapped_phase.content;

  	for (i = 0, i_plus_one = 1, cum = cum_incr = (buf->first_sample_file_pos + sample_index)
  																	* TWO_PI / core_length;
   		 i < last_bin;
       	 i++, i_plus_one++, cum += cum_incr)
    {
        spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
        spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

		aux = c_out[i_plus_one].re * c_out[i_plus_one].re + c_out[i_plus_one].im * c_out[i_plus_one].im;

   		spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   		spec_power_content[i] = aux * POW_NORM_FACTOR;
   		spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;
    }
    /* deals with the problematic last bin separately */
    spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
    spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

    aux = (c_out[i_plus_one].re * c_out[i_plus_one].re + c_out[i_plus_one].im*c_out[i_plus_one].im)/4;

	spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
	spec_power_content[i] = aux * POW_NORM_FACTOR;
	spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;

	/*
   	mag_warp_denoise(spec->power.content, spec->warped_denoised.content,
    													spec->power_noise.content, aud_mod);
	*/
    calc_spectra_extra_data(spec, last_bin + 1);
}

/********************************************************************************************/
/* OLD_NEW_FFT_PH_UPH_MAG_POW_LP_WD_PN                                                      */
/********************************************************************************************/
/* does NOT remove the offset */
void old_new_fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod)
{
	static unsigned char first_call = TRUE;
	static Complex *c_in = NULL, *c_out = NULL;

	unsigned int i, i_plus_one, core_length, total_length, last_bin;
	float *signal;
	double aux, cum, cum_incr, *fft_p_window_content,
			*spec_phase_content, *spec_unwrapped_phase_content, *spec_magnitude_content,
   			*spec_power_content, *spec_log_power_content;

	FFT_properties *fft_p;
	Spectra *spec;

    signal = &(buf->sample_data[channel][sample_index]);
	fft_p = buf->fft_p;
	fft_p_window_content = &(fft_p->window.content[0]);
	spec = &(buf->spec_data[channel][spec_index]);
	core_length = fft_p->window.length;
	total_length = core_length * fft_p->zero_padding_ratio;
	last_bin = fft_p->num_of_bins - 1;

	/*
 	fft_real_signal_to_mag_phase();
 	*/

	if (first_call)
 	{
      c_in = (Complex*) emalloc(total_length * sizeof(Complex));
	  c_out = (Complex*) emalloc(total_length * sizeof(Complex));

   	  for (i = 0; i < core_length; i++)
    	c_in[i].im = 0;

	  for (i = core_length; i < total_length; i++)
	  	c_in[i].re = c_in[i].im = 0;

	  first_call = FALSE;
	}

    if (buf->fft_p->window.apodization_function.type == RECTANGULAR)
    {
    	for (i = 0; i < core_length; i++)
    		c_in[i].re = signal[i];
    }
    else
    {
    	for (i = 0; i < core_length; i++)
    		c_in[i].re = signal[i] * fft_p_window_content[i];
    }

    fft(c_in, c_out, total_length);

    spec_phase_content = spec->phase.content;
    spec_magnitude_content = spec->magnitude.content;
    spec_power_content = spec->power.content;
   	spec_log_power_content = spec->log_power.content;
   	spec_unwrapped_phase_content = spec->unwrapped_phase.content;

  	for (i = 0, i_plus_one = 1, cum = cum_incr = (buf->first_sample_file_pos + sample_index)
  																	* TWO_PI / core_length;
   		 i < last_bin;
       	 i++, i_plus_one++, cum += cum_incr)
    {
        spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
        spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

        aux = c_out[i_plus_one].re * c_out[i_plus_one].re
        						  + c_out[i_plus_one].im * c_out[i_plus_one].im;

   		spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   		spec_power_content[i] = aux * POW_NORM_FACTOR;
   		spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;
    }
    /* deals with the problematic last bin separately */
    spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
    spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

    aux = (c_out[i_plus_one].re * c_out[i_plus_one].re
    						  + c_out[i_plus_one].im * c_out[i_plus_one].im)/4;

	spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
	spec_power_content[i] = aux * POW_NORM_FACTOR;
	spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;

   	mag_warp_denoise(spec->power.content, spec->warped_denoised.content,
    													spec->power_noise.content, aud_mod);
    calc_spectra_extra_data(spec, last_bin + 1);
}

/********************************************************************************************/
/* NEW_FFT_PH_UPH_MAG_POW_LP_WD_PN                                                          */
/********************************************************************************************/
/* DOES remove the offset */
void new_fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod)
{
	static unsigned char first_call = TRUE;
	static Complex *c_in = NULL, *c_out = NULL;

	unsigned int i, i_plus_one, core_length, total_length, last_bin;
	float *signal;
	double aux, cum, cum_incr, *fft_p_window_content,
			*spec_phase_content, *spec_unwrapped_phase_content, *spec_magnitude_content,
   			*spec_power_content, *spec_log_power_content;

	FFT_properties *fft_p;
	Spectra *spec;

    signal = &(buf->sample_data[channel][sample_index]);
	fft_p = buf->fft_p;
	fft_p_window_content = &(fft_p->window.content[0]);
	spec = &(buf->spec_data[channel][spec_index]);
	core_length = fft_p->window.length;
	total_length = core_length * fft_p->zero_padding_ratio;
	last_bin = fft_p->num_of_bins - 1;

	if (first_call)
 	{
      c_in = (Complex*) emalloc(total_length * sizeof(Complex));
	  c_out = (Complex*) emalloc(total_length * sizeof(Complex));

   	  for (i = 0; i < core_length; i++)
    	c_in[i].im = 0;

	  for (i = core_length; i < total_length; i++)
	  	c_in[i].re = c_in[i].im = 0;

	  first_call = FALSE;
	}

	/* LEVEL 1 */

	for (i = 0, aux = 0; i < core_length; i++) {
		aux += signal[i];
	}
	aux /= core_length; /* determine the offset */

	for (i = 0; i < core_length; i++) {
   		c_in[i].re = (signal[i] /*- aux*/) * fft_p_window_content[i]; /* and then removes it */
	}

	/* END OF LEVEL 1 */

	/* LEVEL 2 */
	/* END OF LEVEL 2 */

    fft(c_in, c_out, total_length);

    spec_phase_content = spec->phase.content;
    spec_magnitude_content = spec->magnitude.content;
    spec_power_content = spec->power.content;
   	spec_log_power_content = spec->log_power.content;
   	spec_unwrapped_phase_content = spec->unwrapped_phase.content;

  	for (i = 0, i_plus_one = 1, cum = cum_incr = (buf->first_sample_file_pos + sample_index)
  																	* TWO_PI / core_length;
   		 i < last_bin;
       	 i++, i_plus_one++, cum += cum_incr)
    {
        spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
        spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

        aux = c_out[i_plus_one].re * c_out[i_plus_one].re
        						  + c_out[i_plus_one].im * c_out[i_plus_one].im;

   		spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   		spec_power_content[i] = aux * POW_NORM_FACTOR;
   		spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;
    }
    /* deals with the problematic last bin separately */
    spec_phase_content[i] = atan2(c_out[i_plus_one].im, c_out[i_plus_one].re);
    spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

    aux = (c_out[i_plus_one].re * c_out[i_plus_one].re
        						  + c_out[i_plus_one].im*c_out[i_plus_one].im)/4;

	spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
	spec_power_content[i] = aux * POW_NORM_FACTOR;
	spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;

   	mag_warp_denoise(spec->power.content, spec->warped_denoised.content,
    													spec->power_noise.content, aud_mod);
    calc_spectra_extra_data(spec, last_bin + 1);
}

/********************************************************************************************/
/* PRINT_SPEC                                                                               */
/********************************************************************************************/
void print_spec(Buffer* buf, const unsigned int spec_index, unsigned int offset)
{
  unsigned char c;
  unsigned int b;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate,
    		freq, freq_incr = buf->fft_p->freq_resolution;
  Spectrum *spec_ptr;

  if (buf->stream->num_of_channels == 1)
  {
      spec_ptr = (Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset);
      for (b = 0, freq = freq_incr; b < buf->fft_p->num_of_bins; b++, freq += freq_incr)
      	printf("%g %g %g\n", t, freq, spec_ptr->content[b]);
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
        spec_ptr = (Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset);
      	for (b = 0, freq = freq_incr; b < buf->fft_p->num_of_bins; b++, freq += freq_incr)
          printf("%g %u %g %g\n", t, c + 1, freq, spec_ptr->content[b]);
      }
  }

  /* gnuplot needs this */
  if (DOUBLED_LINEBREAKS)
  	printf("\n\n");
  else
  	printf("\n");
}


/********************************************************************************************/
/* PRINT_SPEC_LOCMAX                                                                        */
/********************************************************************************************/
void print_spec_locmax(Buffer* buf, const unsigned int spec_index, unsigned int offset)
{
  const unsigned int num_of_bins =  buf->fft_p->num_of_bins - 1,
  					loc_max_width = 1+2*buf->fft_p->zero_padding_ratio;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  const double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate,
				freq_res = buf->fft_p->freq_resolution;
  unsigned char c;
  unsigned int k;
  double *spec_content;

  if (buf->stream->num_of_channels == 1)
  {
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k] > spec_content[k - 1] && spec_content[k] > spec_content[k + 1]) {
      		printf("%g %g %g\n", t, freq_res * (1 + k), spec_content[k]);
		}
  } else {
  	for (c = 0; c < buf->stream->num_of_channels; c++)
      {
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
  		for (k = 1; (k < num_of_bins - 1); k++)
      		if (spec_content[k] > spec_content[k - 1] && spec_content[k] > spec_content[k + 1]) {
        		printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + k), spec_content[k]);
      	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}

/********************************************************************************************/
/* PRINT_SPEC_INTERP                                                                        */
/********************************************************************************************/
void print_spec_interp(Buffer* buf, const unsigned int spec_index, unsigned int offset)
{
	/* TODO: implement a generic (method independent) function */
}

/********************************************************************************************/
/* PRINT_SPEC_INTERP_QUINN                                                                  */
/********************************************************************************************/
void print_spec_interp_quinn(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;
  double *spec_content;

  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].magnitude.content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
      	{
      		const double accurate_freq = quinn_2nd_interp(&(buf->spec_data[0][spec_index]), num_of_bins, k),
					accurate_log_power = buf->spec_data[0][spec_index].log_power.content[k] + window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
          	spec_content = buf->spec_data[c][spec_index].magnitude.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = quinn_2nd_interp(&(buf->spec_data[c][spec_index]), num_of_bins, k),
					accurate_log_power = buf->spec_data[c][spec_index].log_power.content[k] + window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}

/********************************************************************************************/
/* PRINT_SPEC_INTERP_GRANDKE                                                                */
/********************************************************************************************/
void print_spec_interp_grandke(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;
  double *spec_content;


  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].magnitude.content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
      	{
      		const double accurate_freq = grandke_interp(buf->spec_data[0][spec_index].magnitude.content, num_of_bins, k),
					accurate_log_power = buf->spec_data[0][spec_index].log_power.content[k] + window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
          	spec_content = buf->spec_data[c][spec_index].magnitude.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = grandke_interp(buf->spec_data[c][spec_index].magnitude.content, num_of_bins, k),
					accurate_log_power = buf->spec_data[c][spec_index].log_power.content[k] + window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);


        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}


/********************************************************************************************/
/* PRINT_SPEC_INTERP_PARABOLIC                                                              */
/********************************************************************************************/
void print_spec_interp_parabolic(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;
  double *spec_content;


  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].magnitude.content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
      	{
      		const double accurate_freq = parabolic_interp(buf->spec_data[0][spec_index].log_power.content, num_of_bins, k),
					accurate_log_power = amplitude_parabolic_interp(buf->spec_data[0][spec_index].log_power.content, num_of_bins, k);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
          	spec_content = buf->spec_data[c][spec_index].magnitude.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = parabolic_interp(buf->spec_data[c][spec_index].log_power.content, num_of_bins, k),
					accurate_log_power = amplitude_parabolic_interp(buf->spec_data[c][spec_index].log_power.content, num_of_bins, k);


        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}


/********************************************************************************************/
/* PRINT_SPEC_INTERP_CENTROID                                                               */
/********************************************************************************************/
void print_spec_interp_centroid(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;
  double *spec_content;


  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].magnitude.content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
      	{
      		const double accurate_freq = centroid_interp(buf->spec_data[0][spec_index].magnitude.content, num_of_bins, k),
					accurate_amplitude = spec_content[k]/hann_win_mag_freq_resp(accurate_freq - (double)k),
					accurate_log_power = 10 * log10(accurate_amplitude * accurate_amplitude);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
          	spec_content = buf->spec_data[c][spec_index].magnitude.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = centroid_interp(buf->spec_data[c][spec_index].magnitude.content, num_of_bins, k),
					accurate_amplitude = spec_content[k]/hann_win_mag_freq_resp(accurate_freq - (double)k),
					accurate_log_power = 10 * log10(accurate_amplitude * accurate_amplitude);


        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}


/********************************************************************************************/
/* PRINT_SPEC_INTERP_DERIVATIVE                                                             */
/********************************************************************************************/
void print_spec_interp_derivative(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;
  double *spec_content;


  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].magnitude.content;
  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
      	{
      		const double accurate_freq = derivative_interp(buf, 0, spec_index, k),
					accurate_amplitude = spec_content[k]/hann_win_mag_freq_resp(accurate_freq - (double)k),
					accurate_log_power = 10 * log10(accurate_amplitude * accurate_amplitude);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
          	spec_content = buf->spec_data[c][spec_index].magnitude.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = derivative_interp(buf, c, spec_index, k),
					accurate_amplitude = spec_content[k]/hann_win_mag_freq_resp(accurate_freq - (double)k),
					accurate_log_power = 10 * log10(accurate_amplitude * accurate_amplitude);


        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}

/********************************************************************************************/
/* PRINT_SPEC_INTERP_REASSIGNMENT                                                           */
/********************************************************************************************/
void print_spec_interp_reassignment(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double *spec_content,
  			t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;

  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].log_power.content;
  	for (k = 1; k < (num_of_bins - 1); k++)

      	if (spec_content[k - 1] <= spec_content[k] && spec_content[k + 1] <= spec_content[k])
      	{
      		const double accurate_freq = reassignment_interp(buf, 0, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
        spec_content = buf->spec_data[c][spec_index].log_power.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = reassignment_interp(buf, c, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}


/********************************************************************************************/
/* GRANDKE_INTERP                                                                           */
/********************************************************************************************/
double grandke_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	double a, d;

	if (bin_to_interp == 0
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp + 1])
		return bin_to_interp;

	if (spec_content[bin_to_interp - 1] > spec_content[bin_to_interp + 1])
	{
		a = spec_content[bin_to_interp]/spec_content[bin_to_interp - 1];
		d =  (2*a - 1)/(a + 1);
		if (fabs(d-1) < 0.52) { /* This comparison is only to avoid letting the interpolation */
			return (bin_to_interp - 1 + d); /* worsen the estimation and 0.52 is empiric. */
		}
	} else {
		a = spec_content[bin_to_interp + 1]/spec_content[bin_to_interp];
		d =  (2*a - 1)/(a + 1);
		if (fabs(d) < 0.52) {
			return (bin_to_interp + d);
		}
	}
	return bin_to_interp;
}

/********************************************************************************************/
/* BARYCENTRIC_INTERP                                                                       */
/********************************************************************************************/
double barycentric_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	double l, c, r;

	if (bin_to_interp == 0
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp + 1])
		return bin_to_interp;

	l = spec_content[bin_to_interp - 1];
	c = spec_content[bin_to_interp];
	r = spec_content[bin_to_interp + 1];

	return ( bin_to_interp + (r - l)/(l + c + r) );
}


/********************************************************************************************/
/* AGREZ_INTERP                                                                             */
/********************************************************************************************/
double agrez_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	double l, c, r, delta;

	if (bin_to_interp == 0
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp + 1])
		return bin_to_interp;

	l = spec_content[bin_to_interp - 1];
	c = spec_content[bin_to_interp];
	r = spec_content[bin_to_interp + 1];

	delta = 2*(r-l)/(l+2*c+r);

	if (fabs(delta) < 0.52) {
		return (bin_to_interp + delta);
	}
	return bin_to_interp;
}

/********************************************************************************************/
/* CHARPENTIER_INTERP                                                                       */
/********************************************************************************************/
/* named after F.J. Charpentier, 86 */
double charpentier_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp)
{
	const double *mag = buf->spec_data[channel][spec_index].magnitude.content,
				*ph = buf->spec_data[channel][spec_index].phase.content,
				w = TWO_PI / buf->fft_p->window.length;
	Complex c, l, r, /* center bin, left and right neighbours */
			a, b;
	double freq_estimate, bin_estimate;

	if (bin_to_interp == 0
		|| mag[bin_to_interp] < mag[bin_to_interp - 1]
		|| bin_to_interp >= (buf->fft_p->num_of_bins - 1)
		|| mag[bin_to_interp] < mag[bin_to_interp + 1])
		return bin_to_interp;

	c_copy(c_polar_to_complex(mag[bin_to_interp - 1], ph[bin_to_interp - 1]), &l);
	c_copy(c_polar_to_complex(mag[bin_to_interp], ph[bin_to_interp]), &c);
	c_copy(c_polar_to_complex(mag[bin_to_interp + 1], ph[bin_to_interp + 1]), &r);

	c_copy(c_mult_by_real(c_add(c, c_mult_by_real(c_add(l, r), -0.5)), 0.5), &a);
	c_copy(c_mult_by_real(c_add(c_mult(c, c_exp(w)), c_mult_by_real(c_add(c_mult(l, c_exp(w)), c_mult(r, c_exp(w))), -0.5)), 0.5), &a);
	/*
	c_copy(c_subtract(c, c_mult_by_real(c_add(c_mult(l, c_exp(w)), c_mult(r, c_exp(w))), 0.5)), &b);
	*/

	/*
	freq_estimate = buf->stream->sample_rate * ((bin_to_interp+1) / (double)buf->fft_p->window.length
													+ INVERSE_OF_TWO_PI * c_phase(c_div(a, b)));
	*/
	freq_estimate = (c_phase(b) - c_phase(a)) * buf->stream->sample_rate;

	bin_estimate = freq_estimate / buf->fft_p->freq_resolution - 1;

	if (1 || fabs(bin_estimate - bin_to_interp) < 0.52) {
		return bin_estimate;
	}
	return bin_to_interp;
}

/********************************************************************************************/
/* PARABOLIC_INTERP                                                                         */
/********************************************************************************************/
double parabolic_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	double l, c, r, delta;

	if (bin_to_interp == 0
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| spec_content[bin_to_interp] < spec_content[bin_to_interp + 1])
		return bin_to_interp;

	l = 20 * log10(spec_content[bin_to_interp - 1]);
	c = 20 * log10(spec_content[bin_to_interp]);
	r = 20 * log10(spec_content[bin_to_interp + 1]);

	delta = (double)1/2*((l - r)/(l - 2*c + r));

	if (fabs(delta) < 0.52) {
		return (bin_to_interp + delta);
	}
	return bin_to_interp;
}

/********************************************************************************************/
/* AMPLITUDE_PARABOLIC_INTERP                                                               */
/********************************************************************************************/
double amplitude_parabolic_interp(double *log_power_spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	double l, c, r, d;

	if (bin_to_interp == 0
		|| log_power_spec_content[bin_to_interp] < log_power_spec_content[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| log_power_spec_content[bin_to_interp] < log_power_spec_content[bin_to_interp + 1])
		return bin_to_interp;

	l = log_power_spec_content[bin_to_interp - 1];
	c = log_power_spec_content[bin_to_interp];
	r = log_power_spec_content[bin_to_interp + 1];

	d = 0.5*(l - r)/(l - 2*c + r);

	return c - (l - r) * d * 0.25;
}


/********************************************************************************************/
/* QUINN_2ND_INTERP                                                                         */
/********************************************************************************************/
double quinn_2nd_interp(Spectra *spectra, const unsigned int num_of_bins,
														const unsigned int bin_to_interp)
{
	const double *mag = spectra->magnitude.content,
					*ph = spectra->phase.content;

	Complex l, c, r;
	double ap, am, dp, dm, d;

	if (bin_to_interp == 0
		|| mag[bin_to_interp] < mag[bin_to_interp - 1]
		|| bin_to_interp >= (num_of_bins - 1)
		|| mag[bin_to_interp] < mag[bin_to_interp + 1])
		return bin_to_interp;

	c_copy(c_polar_to_complex(mag[bin_to_interp - 1], ph[bin_to_interp - 1]), &l);
	c_copy(c_polar_to_complex(mag[bin_to_interp], ph[bin_to_interp]), &c);
	c_copy(c_polar_to_complex(mag[bin_to_interp + 1], ph[bin_to_interp + 1]), &r);

	ap = (r.re*c.re + r.im*c.im)/pow(mag[bin_to_interp], 2);
	dp = ap/(ap - 1);
	am = (l.re*c.re + l.im*c.im)/pow(mag[bin_to_interp], 2);
	dm = am/(1 - am);
	d = 0.5*(dp + dm) + quinn_2nd_interp_polynomial(dp*dp)
					  - quinn_2nd_interp_polynomial(dm*dm);

	if (fabs(d) < 0.52) {
		return bin_to_interp + d;
	}
	return bin_to_interp;
}

/********************************************************************************************/
/* QUINN_2ND_INTERP_POLYNOMIAL                                                              */
/********************************************************************************************/
double quinn_2nd_interp_polynomial(const double x)
{
	return (
			0.25 * log10(3*x*x + 6*x + 1)
			- SQRT_OF_6_DIVIDED_BY_24
				* log10((x + 1 - SQRT_OF_TWO_THIRDS)/(x + 1 + SQRT_OF_TWO_THIRDS))
			);
}

/********************************************************************************************/
/* DERIVATIVE_INTERP                                                                        */
/********************************************************************************************/
/* this implementation is FILTER-ORIENTED, and faster because gain "correction" is applied
 * only for the peak to be interpolated */
/* TODO: fazer rodar direito com zero-padding */
double derivative_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp)
{
	const unsigned int core_len = buf->fft_p->window.length,
					total_len = core_len * buf->fft_p->zero_padding_ratio,
					num_of_bins = total_len/2;
	const float *signal = &(buf->sample_data[channel][buf->fft_p->step * spec_index]);
	const double sample_rate = buf->stream->sample_rate,
			*mag = buf->spec_data[channel][spec_index].magnitude.content;

	static unsigned char first_call = TRUE;
	static unsigned int last_spec_index = UINT_MAX;
	static signed long last_call_first_sample_file_pos = LONG_MAX;
	static double *deriv_mag = NULL;
	static double *deriv_aprox = NULL;

	unsigned int i;
	double accurate_bin;

	if (buf == NULL)
    {
        if (!first_call)
        {
        	if (deriv_aprox) free(deriv_aprox);
        	if (deriv_mag) free(deriv_mag);

        	first_call = TRUE;
        }

        return -1;
    }

	if (first_call)
	{
		deriv_mag = (double*)emalloc(num_of_bins * sizeof(double));
		deriv_aprox = (double*)emalloc(core_len * sizeof(double));

		first_call = FALSE;
	}

	if (spec_index != last_spec_index
		|| last_call_first_sample_file_pos != buf->first_sample_file_pos)
	{
		for (i = 0; i < core_len; i++) {
			deriv_aprox[i] = (signal[i] - signal[i - 1]) * sample_rate;
		}

		fft_real_signal_to_mag(deriv_aprox, buf->fft_p, deriv_mag);

		last_spec_index = spec_index;
		last_call_first_sample_file_pos = buf->first_sample_file_pos;
	}

	accurate_bin  = (sample_rate*INVERSE_OF_PI * asin(deriv_mag[bin_to_interp]/(2*sample_rate*mag[bin_to_interp]))) / buf->fft_p->freq_resolution - 1;

	if (fabs(accurate_bin - bin_to_interp) < 0.52) /* valor empírico, mas parece que não é sequer necessário restringir
													  exceto em casos em que os pressupostos não valem (coincidência, por ex) */
		return accurate_bin;

	return bin_to_interp;
}

/********************************************************************************************/
/* SLOWER_DERIVATIVE_INTERP                                                                 */
/********************************************************************************************/
/* this implementation is DERIVATIVE-ORIENTED, and slower because gain correction
 * is applied to the whole spectrum */
double slower_derivative_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp)
{
	const unsigned int core_len = buf->fft_p->window.length,
					total_len = core_len * buf->fft_p->zero_padding_ratio,
					num_of_bins = total_len/2;
	const float *signal = &(buf->sample_data[channel][buf->fft_p->step * spec_index]);
	const double sample_rate = buf->stream->sample_rate,
			*mag = buf->spec_data[channel][spec_index].magnitude.content;

	static unsigned char first_call = TRUE;
	static unsigned int last_spec_index = UINT_MAX;
	static double *deriv_mag = NULL;

	unsigned int i;
	double deriv_aprox[core_len], accurate_freq, accurate_bin, ang_freq, gain_factor;

	if (first_call)
	{
		deriv_mag = (double*)emalloc(num_of_bins * sizeof(double));
		first_call = FALSE;
	}

	if (spec_index != last_spec_index)
	{
		for (i = 0; i < core_len; i++) {
			deriv_aprox[i] = (signal[i] - signal[i - 1]) * sample_rate;
		}

		fft_real_signal_to_mag(deriv_aprox, buf->fft_p, deriv_mag);

		for (i = 0; i < num_of_bins; i++) {
			ang_freq = (double)(i + 1) / num_of_bins * PI;
			gain_factor = ang_freq / (2 * sin(ang_freq/2));
			deriv_mag[i] *= gain_factor;
		}

		last_spec_index = spec_index;
	}

	accurate_freq = (deriv_mag[bin_to_interp]/mag[bin_to_interp]) * INVERSE_OF_TWO_PI;
	accurate_bin = accurate_freq / buf->fft_p->freq_resolution - 1;

	if (fabs(accurate_bin - bin_to_interp) < 0.5) /* menor ou menor-ou-igual, eis a questao? */
		return accurate_bin;

	return bin_to_interp;
}

/********************************************************************************************/
/* REASSIGNMENT_INTERP                                                                      */
/********************************************************************************************/
double reassignment_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp)
{
	const unsigned int core_len = buf->fft_p->window.length,
					total_len = core_len * buf->fft_p->zero_padding_ratio;
	const float *signal = &(buf->sample_data[channel][buf->fft_p->step * spec_index]);

	static unsigned char first_call = TRUE;
	static unsigned int last_spec_index = UINT_MAX;
	static signed long last_call_buf_first_sample_file_pos = LONG_MAX;
	static Complex *c_dh_win_spec = NULL;
	static double *windowed_signal = NULL;

	unsigned int i;
	double accurate_bin;

	if (buf == NULL)
    {
        if (!first_call)
        {
        	if (c_dh_win_spec) free(c_dh_win_spec);
        	if (windowed_signal) free(windowed_signal);

        	first_call = TRUE;
        }

        return 0;
    }

	if (first_call == TRUE)
	{
		c_dh_win_spec = (Complex*)emalloc(total_len/2 * sizeof(Complex));
		windowed_signal = (double*)emalloc(total_len * sizeof(double));

		for (i = core_len; i < total_len; i++) {
			windowed_signal[i] = 0;
		}


		first_call = FALSE;
	}

	if (spec_index != last_spec_index
		|| last_call_buf_first_sample_file_pos != buf->first_sample_file_pos)
	{
		for (i = 0; i < core_len; i++) {
			windowed_signal[i] = (signal[i] * buf->fft_p->window.derivative[i]);
		}
		fft_real_signal(c_dh_win_spec, windowed_signal, total_len);

		last_spec_index = spec_index;
		last_call_buf_first_sample_file_pos = buf->first_sample_file_pos;
	}

	accurate_bin = bin_to_interp + (total_len/TWO_PI) * c_div(
			c_mult_by_real(c_dh_win_spec[bin_to_interp], MAG_NORM_FACTOR),
			c_polar_to_complex(buf->spec_data[channel][spec_index].magnitude.content[bin_to_interp],
								buf->spec_data[channel][spec_index].phase.content[bin_to_interp]) ).im;

	if (fabs(accurate_bin - bin_to_interp) < 0.52) /* valor empírico, copiado do derivative_interp */
		return accurate_bin;

	return bin_to_interp;
}


/********************************************************************************************/
/* CENTROID_INTERP                                                                          */
/********************************************************************************************/
double centroid_interp(double *v, const unsigned int n, const unsigned int center)
{
	unsigned int const max_neightbours = 1;
    int bound_dist, neightbours_to_be_used, i;
    double abs_total, total, average;
    Window w;

    if (center == 0 || center == (n - 1) )
        return center;

    bound_dist = (center < n/2 ? center : n - center - 1);

    for (i = 0; i <= bound_dist; i++)
    {
        if (v[center - i - 1] >= v[center - i]) break;
        if (v[center + i] <= v[center + i + 1]) break;
    }

    neightbours_to_be_used = (i - 1);

    if (neightbours_to_be_used > max_neightbours)
        neightbours_to_be_used = max_neightbours;

    if (neightbours_to_be_used == 0)
        return center;

    w.length = (2 * neightbours_to_be_used + 1);
    w.content = (double *) malloc (sizeof(double) * w.length);

    if (w.length == 3)
    {
        w.content[0] = 3.9;
        w.content[1] = 5.0;
        w.content[2] = w.content[0];
    }

    for (abs_total = w.content[neightbours_to_be_used] * v[center], total = 0.0, i = 1; i <= neightbours_to_be_used; i++)
    {
        total += ( w.content[neightbours_to_be_used + i] * v[center + i]
        			- w.content[neightbours_to_be_used - i] * v[center - i] );
        abs_total += (v[center - i] + v[center + i]);
    }

    average = (double)neightbours_to_be_used * total / abs_total;

    return (center + average);
}


/********************************************************************************************/
/* PRINT_UNPRED                                                                             */
/********************************************************************************************/
void print_unpred(Buffer* buf, unsigned int index)
{
    double t = ((double)buf->first_sample_file_pos + buf->unpredictability_to_sample_offset
    			+ index * buf->fft_p->step)	/ buf->stream->sample_rate;

	printf("%g %g\n", t, buf->unpredictability_data[index]);
}

/********************************************************************************************/
/* PRINT_THRESHOLD                                                                          */
/********************************************************************************************/
void print_threshold(Buffer* buf, unsigned int index)
{
    double t = ((double)buf->first_sample_file_pos + buf->threshold_to_sample_offset
    			+ index * buf->fft_p->step) / buf->stream->sample_rate;

	printf("%g %g\n", t, buf->threshold_data[index]);
}

/********************************************************************************************/
/* INIT_AUDITORY_MODEL                                                                      */
/********************************************************************************************/
void init_auditory_model(Auditory_Model **auditory_model, Buffer* buf)
{
    unsigned char i, num_of_cb;
    unsigned int j, num_of_bins = buf->fft_p->num_of_bins;
    double freq_res = buf->fft_p->freq_resolution, w, f0_min = max(MIN_ABSOLUTE_F0,
    					MIN_F0_TO_FREQ_RES_RATIO * freq_res * buf->fft_p->zero_padding_ratio),
         	temp_band_center, temp_band_width, temp_band_response_total_sum;

    Auditory_Model *aud_mod;

    if (buf->stream->sample_rate < MAX_ABSOLUTE_F0 * 2)
    {
        fprintf(stderr, "Error in init_auditory_model: the sample rate is too low.\n");
        exit(ERROR); /* TODO: change void to signed int and avoid using exit() */
    }

    if (!*auditory_model) *auditory_model = (Auditory_Model*) emalloc(sizeof(Auditory_Model));

    aud_mod = *auditory_model;

    temp_band_width = max(log2((f0_min + 100) / f0_min), (double)2/3);
    temp_band_center = f0_min * pow(2, temp_band_width/2);
    num_of_cb = 0;
    do  /* TODO: ensure that the case in which no C.B. can be created is well handled */
    {
        num_of_cb++;
        temp_band_width = max(log2((temp_band_center + 100)
     									/ temp_band_center), (double)2/3);
      	temp_band_center = temp_band_center * pow(2, temp_band_width/2);
    } while (num_of_cb < MAX_CRITICAL_BANDS
    		 && temp_band_center * pow(2, temp_band_width/2) < MAX_ABSOLUTE_F0);

    aud_mod->fft_p = buf->fft_p;
    aud_mod->num_of_bands = num_of_cb;
    aud_mod->band = (Critical_Band*) emalloc(num_of_cb * sizeof(Critical_Band));

    aud_mod->band[0].response = (double*) ecalloc(num_of_bins * sizeof(double));
    aud_mod->band[0].width = max(log2((f0_min + 100) / f0_min), (double)2/3);
    aud_mod->band[0].center = f0_min * pow(2, aud_mod->band[0].width/2);

   	for (i = 1; i < num_of_cb; i++)
    {
    	aud_mod->band[i].response = (double*) ecalloc(num_of_bins * sizeof(double));
    	aud_mod->band[i].width = max(log2((aud_mod->band[i - 1].center + 100)
     									/ aud_mod->band[i - 1].center), (double)2/3);
      	aud_mod->band[i].center = aud_mod->band[i - 1].center
       													* pow(2, aud_mod->band[i].width/2);
   	}

    for (i = 0; i < num_of_cb; i++)
    {
        w = pow(2, aud_mod->band[i].width / 2);
        aud_mod->band[i].first_bin = ceil(aud_mod->band[i].center / w / freq_res) - 1;
        aud_mod->band[i].num_of_bins = (floor(aud_mod->band[i].center * w / freq_res) - 1)
        								- aud_mod->band[i].first_bin + 1;

        for (j = aud_mod->band[i].first_bin, temp_band_response_total_sum = 0.0;
        	 j < aud_mod->band[i].first_bin + aud_mod->band[i].num_of_bins;
          	 j++)
        {
            if (freq_res * (j + 1) <= aud_mod->band[i].center) {

                aud_mod->band[i].response[j] = (log((j + 1) * freq_res)
                								- log(aud_mod->band[i].center))/log(w) + 1;
            }
            else {

                aud_mod->band[i].response[j] = -1 * (log((j + 1) * freq_res)
                								- log(aud_mod->band[i].center))/log(w) + 1;
            }

            temp_band_response_total_sum += aud_mod->band[i].response[j];

       	}

       aud_mod->band[i].response_total_sum = temp_band_response_total_sum;
    }

    return;
}

/********************************************************************************************/
/* PRINT_AUDITORY_MODEL                                                                     */
/********************************************************************************************/
void print_auditory_model(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;
    unsigned int k;
    double freq_res = aud_mod->fft_p->freq_resolution;

    for (b = 0; b < num_of_bands; b++)
    {
    	for (k = 0; k <  aud_mod->fft_p->num_of_bins; k++)
          printf("%g %u %g\n", (k + 1) * freq_res, b, aud_mod->band[b].response[k]);

		printf("\n\n");
    }

    return;
}

/********************************************************************************************/
/* PRINT_AUDITORY_MODEL_SUMMARY                                                             */
/********************************************************************************************/
void print_auditory_model_summary(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;
    unsigned int k, last_k, num_of_bins = aud_mod->fft_p->num_of_bins;
    double freq_res = aud_mod->fft_p->freq_resolution;

    for (b = 0; b < num_of_bands; b++)
    {
    	for (k = max((aud_mod->band[b].first_bin - 1), 0),
     		 	last_k = min((aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins),
        																		num_of_bins);
     		 k <= last_k;
        	 k++)
      	{
          printf("%g %u %g\n", (k + 1) * freq_res, b, aud_mod->band[b].response[k]);
        }

		printf("\n");
    }

    return;
}

/********************************************************************************************/
/* PRINT_AUDITORY_MODEL_PER_BIN_SUM                                                         */
/********************************************************************************************/
void print_auditory_model_per_bin_sum(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;
    unsigned int k, num_of_bins = aud_mod->fft_p->num_of_bins;
    double freq_res = aud_mod->fft_p->freq_resolution, sum;

    for (k = 0; k < num_of_bins; k++)
    {
        sum = 0;
    	for (b = 0; b < num_of_bands; b++)
    		sum += aud_mod->band[b].response[k];

    	printf("%g %g\n", (k + 1) * freq_res, sum);
  	}

    return;
}

/********************************************************************************************/
/* PRINT_AUDITORY_MODEL_PER_BAND_SUM                                                        */
/********************************************************************************************/
void print_auditory_model_per_band_sum(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;

    for (b = 0; b < num_of_bands; b++)
    	printf("%u %g\n", (b + 1), aud_mod->band[b].response_total_sum);

    return;
}

/********************************************************************************************/
/* PRINT_AUDITORY_MODEL_PER_BAND_AVERAGE                                                    */
/********************************************************************************************/
void print_auditory_model_per_band_average(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;

    for (b = 0; b < num_of_bands; b++)
    	printf("%u %g\n", (b + 1), aud_mod->band[b].response_total_sum
     														/ aud_mod->band[b].num_of_bins);

    return;
}


/********************************************************************************************/
/* PRINT_AUDITORY_MODEL_LIMITS                                                              */
/********************************************************************************************/
void print_auditory_model_limits(Auditory_Model* aud_mod)
{
    unsigned char b, num_of_bands = aud_mod->num_of_bands;
    double freq_res = aud_mod->fft_p->freq_resolution;

    for (b = 0; b < num_of_bands; b++)
    {
    	printf("%u %g %g\n", b, aud_mod->band[b].first_bin * freq_res,
				(aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins - 1) * freq_res);
  	}

    return;
}


/* Main F0 methods to be compared with the proposed:
 * Time domain:
 *  YIN
 *  AC (Eq.1)
 *  AMDF
 *  WAC: AC(n)/(AMDF(n)+k), typically k=1
 *  MWAC (THESIS.PDF, Eq.3.7, pg 27)
 * Frequency domain:
 *  CEPSTRUM
 *  HPS

 *  Partial domain:
 *   Piszczalski & Galler
 * 	 DWS (Goldstein)
 *   Hermes (Subharmonic summation)
 *   Brown & Puckette (High-resolution ... based on phase changes)
 *
 * */

/********************************************************************************************/
/* YIN_F0_ESTIMATE                                                                          */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void yin_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	const double THRESHOLD = 0.1;
	const unsigned int quar_fr_len = frame_len/4, half_fr_len = frame_len/2,
						max_lag = min(half_fr_len, MAX_ALLOWED_LAG);
	unsigned int lag, i, i_max;
	double *df, *cmndf, cum;

	*period_estimate = 0;
	*estimate_quality = 1;

	df = (double*)emalloc((max_lag+1)*sizeof(double));

	df[0] = 0;
	for (lag = 1; lag <= max_lag; lag++) {
		double d = 0;
		for (i = quar_fr_len - (lag>>1), i_max = i + half_fr_len - 1;
			 i <= i_max;
			 i++)
		{
			const double p = s[i]-s[i+lag];
			d += p*p;
		}
		df[lag] = d;
	}

	cmndf = (double*)emalloc((max_lag+1)*sizeof(double));

	cmndf[0] = 1;
	for (lag = 1, cum = 0; lag <= max_lag; lag++) {
		cum += df[lag];
		cmndf[lag] = df[lag]*lag/cum;
	}
	free(df);

	for (lag = MIN_ALLOWED_LAG; lag < max_lag; lag++) {
		double l, c, r, period_delta, parab_cmndf;

		if (cmndf[lag] < cmndf[lag-1] && cmndf[lag] < cmndf[lag+1]) {
			l = cmndf[lag-1];
			c = cmndf[lag];
			r = cmndf[lag+1];
			period_delta = 0.5*(l - r)/(l - 2*c + r);
			parab_cmndf = c - (l - r) * period_delta * 0.25;

			if (parab_cmndf < *estimate_quality) {
				l = df[lag-1];
				c = df[lag];
				r = df[lag+1];
				period_delta = 0.5*(l - r)/(l - 2*c + r);

				*period_estimate = lag + period_delta;
				*estimate_quality = parab_cmndf;
			}

			if (parab_cmndf < THRESHOLD) break;
		}
	}
	free(cmndf);
}


/********************************************************************************************/
/* YAN_F0_ESTIMATE                                                                          */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void yan_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	const double THRESHOLD = 0.1;
	const unsigned int max_lag = min(frame_len/2, MAX_ALLOWED_LAG);
	unsigned int lag, i;
	double *df, *cmndf, cum;

	*period_estimate = 0;
	*estimate_quality = 1;

	df = (double*)emalloc((max_lag+1)*sizeof(double));

	df[0] = 0;
	for (lag = 1; lag <= max_lag; lag++) {
		double d = 0;
		for (i = 0; i < max_lag; i++) {
			const double p = s[i]-s[i+lag];
			d += (p >= 0 ? p : -p);
		}
		df[lag] = d;
	}

	cmndf = (double*)emalloc((max_lag+1)*sizeof(double));

	cmndf[0] = 1;
	for (lag = 1, cum = 0; lag <= max_lag; lag++) {
		cum += df[lag];
		cmndf[lag] = df[lag]*lag/cum;
	}
	free(df);

	for (lag = MIN_ALLOWED_LAG; lag < max_lag; lag++) {
		double l, c, r, period_delta, parab_cmndf;

		if (cmndf[lag] < cmndf[lag-1] && cmndf[lag] < cmndf[lag+1]) {
			l = cmndf[lag-1];
			c = cmndf[lag];
			r = cmndf[lag+1];
			period_delta = 0.5*(l - r)/(l - 2*c + r);
			parab_cmndf = c - (l - r) * period_delta * 0.25;

			if (parab_cmndf < *estimate_quality) {
				l = df[lag-1];
				c = df[lag];
				r = df[lag+1];
				period_delta = 0.5*(l - r)/(l - 2*c + r);

				*period_estimate = lag + period_delta;
				*estimate_quality = parab_cmndf;
			}

			if (parab_cmndf < THRESHOLD) break;
		}
	}
	free(cmndf);
}


/********************************************************************************************/
/* GET_MAX_THR_INTERP                                                                           */
/********************************************************************************************/
void get_max_thr_interp(const double df[], const unsigned int max_lag,
						const double threshold,
						double *period_estimate, double *estimate_quality)
{
	unsigned int lag, max_ac_lag = 0;
	for (lag = MIN_ALLOWED_LAG; lag < max_lag; lag++) {
		if ((max_ac_lag == 0 || df[lag] > df[max_ac_lag])
			&& df[lag] > df[lag-1] && df[lag] > df[lag+1]) {
			max_ac_lag = lag;
		}
		if (max_ac_lag != 0 && df[max_ac_lag] > threshold*df[0]) break;
	}

	if (max_ac_lag == 0) {
		max_ac_lag = MIN_ALLOWED_LAG;
		for (lag = MIN_ALLOWED_LAG+1; lag < max_lag; lag++) {
			if (df[lag] > df[max_ac_lag]) {
				max_ac_lag = lag;
			}
		}
		*period_estimate = max_ac_lag;
		*estimate_quality = df[max_ac_lag]/df[0];
	} else {
		double l, c, r, period_delta, parab_df;

		l = df[max_ac_lag-1];
		c = df[max_ac_lag];
		r = df[max_ac_lag+1];
		if (l < c && c > r) {
			period_delta = 0.5*(l - r)/(l - 2*c + r);
			parab_df = c - (l - r) * period_delta * 0.25;

			*period_estimate = max_ac_lag + period_delta;
			*estimate_quality = parab_df/df[0];
		} else {
			*period_estimate = max_ac_lag;
			*estimate_quality = df[max_ac_lag]/df[0];
		}
	}
	*estimate_quality = min(*estimate_quality, 1); /* avoids unpleasant values > 1 */
}

/********************************************************************************************/
/* GET_MIN_THR_INTERP                                                                           */
/********************************************************************************************/
void get_min_thr_interp(const double df[], const unsigned int max_lag,
						const double threshold,
						double *period_estimate, double *estimate_quality)
{
	unsigned int lag, max_ac_lag = 0;
	for (lag = MIN_ALLOWED_LAG; lag < max_lag; lag++) {
		if ((max_ac_lag == 0 || df[lag] < df[max_ac_lag])
			&& df[lag] < df[lag-1] && df[lag] < df[lag+1]) {
			max_ac_lag = lag;
		}
		if (max_ac_lag != 0 && df[max_ac_lag] < threshold) break;
	}

	if (max_ac_lag == 0) {
		max_ac_lag = MIN_ALLOWED_LAG;
		for (lag = MIN_ALLOWED_LAG+1; lag < max_lag; lag++) {
			if (df[lag] < df[max_ac_lag]) {
				max_ac_lag = lag;
			}
		}
		*period_estimate = max_ac_lag;
		*estimate_quality = df[max_ac_lag];
	}
	double l, c, r, period_delta, parab_df;

	l = df[max_ac_lag-1];
	c = df[max_ac_lag];
	r = df[max_ac_lag+1];
	if (l > c && c < r) {
		period_delta = 0.5*(l - r)/(l - 2*c + r);
		parab_df = c - (l - r) * period_delta * 0.25;

		*period_estimate = max_ac_lag + period_delta;
		*estimate_quality = parab_df/df[0];
	}
}

void remove_offset(float s[], const unsigned int n)
{
	unsigned int i;
	double mean = 0;

	for (i = 0; i < n; i++) {
		mean += s[i];
	}
	mean /= n;
	for (i = 0; i < n; i++) {
		s[i] -= mean;
	}
}

/********************************************************************************************/
/* AMDF_F0_ESTIMATE                                                                          */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void amdf_f0_estimate(const float sig[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	const double THRESHOLD = 1.0;
	const unsigned int max_lag = min(frame_len/2, MAX_ALLOWED_LAG);
	unsigned int lag, i;
	double *df;
	float *s;

	*period_estimate = 0;
	*estimate_quality = 1;

	df = (double*)emalloc((max_lag+1)*sizeof(double));

	s = (float*)emalloc(frame_len*sizeof(double));
	memmove(s, sig, frame_len*sizeof(float));

	remove_offset(s, frame_len);

	for (lag = 0; lag <= max_lag; lag++) {
		double d = 0;
		for (i = 0; i < max_lag; i++) {
			const double p = s[i]-s[i+lag];
			d += (p >= 0 ? p : -p);
		}
		df[lag] = d/max_lag;
	}
	get_min_thr_interp(df, max_lag, THRESHOLD, period_estimate, estimate_quality);

	free(df);
	free(s);
}


/********************************************************************************************/
/* AC_F0_ESTIMATE                                                                           */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void ac_f0_estimate(const float sig[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	const double THRESHOLD = 0.9;
	const unsigned int max_lag = min(frame_len/2, MAX_ALLOWED_LAG),
					half_fr_len = frame_len/2, quar_fr_len = frame_len/4;
	unsigned int lag, i, i_max;
	double *df;
	float *s;

	*period_estimate = 0;
	*estimate_quality = 1;

	df = (double*)emalloc((max_lag+1)*sizeof(double));

	s = (float*)emalloc(frame_len*sizeof(double));
	memmove(s, sig, frame_len*sizeof(float));

	remove_offset(s, frame_len);

	for (lag = 0; lag <= max_lag; lag++) {
		double d = 0;
		for (i = quar_fr_len - (lag>>1), i_max = i + half_fr_len - 1;
			 i <= i_max;
			 i++)
		{
			d += s[i]*s[i+lag];
		}
		df[lag] = d/max_lag;
	}
	get_max_thr_interp(df, max_lag, THRESHOLD, period_estimate, estimate_quality);


	free(df);
	free(s);
}

/********************************************************************************************/
/* WAC_F0_ESTIMATE                                                                          */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void wac_f0_estimate(const float sig[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	const double k = 1, /* for stabilization */
				THRESHOLD = 1;
	const unsigned int max_lag = min(frame_len/2, MAX_ALLOWED_LAG);
	unsigned int lag, i;
	double *df;
	float *s;

	*period_estimate = 0;
	*estimate_quality = 0;

	df = (double*)emalloc((max_lag+1)*sizeof(double));

	s = (float*)emalloc(frame_len*sizeof(double));
	memmove(s, sig, frame_len*sizeof(float));

	for (lag = 0; lag <= max_lag; lag++) {
		double ac = 0, amdf = 0;
		for (i = 0; i < max_lag; i++) {
			const double p = s[i]-s[i+lag];
			amdf += (p >= 0 ? p : -p);
			ac += s[i]*s[i+lag];
		}
		df[lag] = ac/(amdf+k);
	}
	get_max_thr_interp(df, max_lag, THRESHOLD, period_estimate, estimate_quality);

	free(df);
	free(s);
}

/********************************************************************************************/
/* CEPSTRUM_F0_ESTIMATE                                                                           */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void cepstrum_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	double *pow;
	const unsigned int ZPF = 1, n = frame_len*ZPF, num_of_bins = n/2+1;


	pow = (double*)emalloc(num_of_bins*sizeof(double));

	reference_fft_real_signal_power_phase(s, frame_len, HAMMING, n-frame_len, pow, NULL);


	float *lp = (float*)emalloc(n*sizeof(float));
	const unsigned int center = num_of_bins-2;
	unsigned int i;
	for (i = 0; i < num_of_bins-1; i++) {
		lp[center+i] = lp[center-i] = log10(pow[i]);
	}
	lp[center+num_of_bins-1] = log10(pow[num_of_bins-1]);
	remove_offset(lp, n);

	reference_fft_real_signal_power_phase(lp, frame_len, RECTANGULAR, n-frame_len, pow, NULL);

	unsigned int max_i = 0;
	const unsigned int max_lag = min(num_of_bins-1, MAX_ALLOWED_LAG);
	for (i = MIN_ALLOWED_LAG; i < max_lag; i++) {
		if (pow[i] > pow[i-1] && pow[i] > pow[i+1]
			&& (max_i == 0 || pow[i] > pow[max_i])) {
			max_i = i;
		}
	}

	*period_estimate = max_i;
	*estimate_quality = 1;

	free(pow);
	free(lp);
}

/********************************************************************************************/
/* ACLOS_F0_ESTIMATE                                                                           */
/********************************************************************************************/
/* requirements: frame_len must be even */
/* note that not all frame is used for all lags which causes delay in f0 estimation */
void aclos_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality)
{
	if (frame_len%2) {
		*period_estimate = DBL_MIN;
		*estimate_quality = DBL_MIN;
		return;
	}

	double *pow;
	float *lp;
	const unsigned int ZPF = 1, n = ZPF*frame_len, num_of_bins = n/2;
	unsigned int i;

	pow = (double*)emalloc(num_of_bins*sizeof(double));
	lp = (float*)emalloc(num_of_bins*sizeof(float));

	reference_fft_real_signal_power_phase(s, frame_len, HANN, n-frame_len, pow, NULL);

	for (i = 1; i < num_of_bins; i++) {
		lp[i-1] = log10(pow[i]);
	}

	const double THRESHOLD = 0.9;
	double *df;
	unsigned int lag, max_lag = num_of_bins/2;

	df = (double*)emalloc(num_of_bins*sizeof(double));

	remove_offset(lp, num_of_bins);

	for (lag = 0; lag <= max_lag; lag++) {
		double d = 0;
		for (i = 0; i < max_lag; i++) {
			d += lp[i]*lp[i+lag];
		}
		df[lag] = d/max_lag;
	}
	get_max_thr_interp(df, max_lag, THRESHOLD, period_estimate, estimate_quality);


	free(df);

	*period_estimate = frame_len / (*period_estimate * ZPF);
	/*
	ac_f0_estimate(lp+1, num_of_bins-1, period_estimate, estimate_quality);
	ac_f0_estimate(s, frame_len, period_estimate, estimate_quality);
	*/

	free(pow);
	free(lp);
}


/********************************************************************************************/
/* WIDEBAND_PITCH_ESTIMATE                                                                  */
/********************************************************************************************/
void klapuri_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned int f0_max_bin = aud_mod->band[aud_mod->num_of_bands - 1].first_bin
    								+ aud_mod->band[aud_mod->num_of_bands - 1].num_of_bins - 1;
    double **salience, aux_salience, aux_sum,
    		f0_min = aud_mod->band[0].first_bin * aud_mod->fft_p->freq_resolution;
    unsigned char b;
    unsigned int h, j, k, k0, k1, m, m0, m1, n, n0, n1, lb, max_offset, aux_index;

    salience = (double**) emalloc(aud_mod->num_of_bands * sizeof(double*));

    for (b = 0; b < aud_mod->num_of_bands; b++);
   	{
        salience[b] = (double*) emalloc((f0_max_bin + 1)* sizeof(double));
    	for (k = 0; k <= f0_max_bin; k++)
      		salience[b][k] = 0;
    }

    for (b = 0; b < aud_mod->num_of_bands; b++)
    {
        printf("bla 8: %d\n", b);
        n0 = floor(f0_min / aud_mod->fft_p->freq_resolution);
        n1 = aud_mod->band[b].num_of_bins - 1;
        lb = aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins - 1;
        for (n = n0; n <= n1; n++)
        {
            printf("bla 8b: %d\n", b);
            m0 = floor(ceil((double)aud_mod->band[b].first_bin / n) * n + 0.5)
            													- aud_mod->band[b].first_bin;
            max_offset = lb * aud_mod->fft_p->freq_resolution
            				* (sqrt(1 + 0.01 * (pow((double)lb/n, 2) - 1)) - 1);
            m1 = m0 + max_offset;
            if (m1 > m0 + n - 1)
            {
                m0 = 0;
                m1 = n - 1;
            }
            salience[b][n] = 0; /* unnecessary, since the matrix was initialized above */
            for (m = m0; m <= m1; m++)
            {
                printf("bla 8c: %d\n", b);
                j = floor((aud_mod->band[b].num_of_bins - m - 1)/n) + 1;
                for (aux_index = 0, aux_sum = 0; aux_index <= j - 1; j++)
                {
                    printf("bla 8d: %d\n", aud_mod->band[b].first_bin + m + n * aux_index);
                	aux_sum += spec->content[aud_mod->band[b].first_bin + m + n * aux_index];
               	}
                aux_salience = (0.75 / j + 0.25) * aux_sum;
                if (aux_salience > salience[b][n])
                	salience[b][n] = aux_salience;
            }
            printf("bla 8e: %d\n", b);
        }

        h = 1;
        k0 = floor( (aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins) / (h+1) );
        if (k0 < aud_mod->band[b].first_bin)
        	k0 = aud_mod->band[b].first_bin;
        k1 = aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins - 1;
        while (k0 <= k1)
        {
            for (k = k0; k <= k1; k++)
            {
                n = floor(k/h + 0.5);
                if (salience[b][n] < spec->content[k])
                	salience[b][n] = spec->content[k];
            }

            h = h + 1;
            k0 = ceil( (aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins) * h
            																		/ (h+1) );
            if (k0 < aud_mod->band[b].first_bin)
            	k0 = aud_mod->band[b].first_bin;
            k1 = floor( (aud_mod->band[b].first_bin - 1) * h / (h - 1) );
            if (k1 > aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins)
            	k1 = aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins;
        }
    }

    printf("bla 9\n");

    for (b = 0; b < aud_mod->num_of_bands; b++)
    	for (n = 0; n <= f0_max_bin; n++)
    		printf("%u %g %g\n", b, (n + 1) * aud_mod->fft_p->freq_resolution,
																			salience[b][n]);

}

/********************************************************************************************/
/* MAX_INDEX_PITCH_ESTIMATE                                                                 */
/********************************************************************************************/
/* currently doest NOT use any kind interpolation */
void max_index_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
   const double strongest_freq = (spec->max_index + 1) * aud_mod->fft_p->freq_resolution,
   					note0_fund_freq = pr->note[0].fund_freq;
   unsigned char max_note;

   if (strongest_freq >= note0_fund_freq)
   {
       max_note = 12 * (int)(12 * log2(strongest_freq/note0_fund_freq) + 0.5);

       if (max_note < pr->num_of_notes)
       {
       		pr->note[max_note].last_evidence_time = curr_time;
			if (OPERATION_MODE == PRINT_F0_ESTIMATE)
				printf("%g %u 1\n", curr_time, max_note);
       }
   }

   return;
}

/********************************************************************************************/
/* HPS_PITCH_ESTIMATE                                                                       */
/********************************************************************************************/
void hps_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned char note, h, num_of_harm = 4, max_note_index = 0, harm[4] = {1, 2, 3, 4},
    	num_of_notes = pr->num_of_notes, first_note_midi_number = pr->note[0].midi_number;;
    double hps, max_val = (-1) * DBL_MAX, fund_freq,
    		freq_res = aud_mod->fft_p->freq_resolution;

    for (note = 0; note < num_of_notes; note++)
    {
    	for (h = 0, hps = 1, fund_freq = pr->note[note].fund_freq;
     		 h < num_of_harm;
        	 h++)
      	{
    		hps *= spec->content[(int)(((fund_freq * harm[h]) / freq_res) + 0.5) - 1];
        }

        if (OPERATION_MODE == PRINT_F0_ESTIMATE)
			printf("%g %u %g\n", curr_time, (note + first_note_midi_number), hps);

        else if (hps > max_val)
     	{
     		max_val = hps;
          	max_note_index = note;
    	}
    }
    /* gnuplots needs this */
    if (OPERATION_MODE == PRINT_F0_ESTIMATE) printf("\n");
    else pr->note[max_note_index].last_evidence_time = curr_time;

    return;
}

/********************************************************************************************/
/* HSC_PITCH_ESTIMATE                                                                       */
/********************************************************************************************/
void hsc_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned char note, h, num_of_harm = 7, max_note_index = 0,
    				num_of_notes = pr->num_of_notes;
    unsigned char *harm;
    unsigned int harm_bin, num_of_bins = aud_mod->fft_p->num_of_bins;
    double hsc, max_val = (-1) * DBL_MAX, fund_freq,
    		freq_res = aud_mod->fft_p->freq_resolution;

    harm = (unsigned char*)emalloc(num_of_harm * sizeof(unsigned char));

    for (h = 0; h < num_of_harm; h++)
    	harm[h] = (h + 1);

    for (note = 0; note < num_of_notes; note++)
    {
        fund_freq = pr->note[note].fund_freq;

        if (fund_freq < freq_res)
        	continue;

    	for (h = 0, hsc = 0; h < num_of_harm; h++)
      	{
           harm_bin = (int)(((fund_freq * harm[h]) / freq_res) + 0.5) - 1;

           if (harm_bin < num_of_bins)
           		hsc += (spec->content[harm_bin] / harm[h]);
           else break;
        }

     	printf("%g %u %g\n", curr_time, (note + 1), hsc);

     	if (hsc > max_val)
     	{
     		max_val = hsc;
          	max_note_index = note;
    	}
    }
    /* gnuplots needs this */
    printf("\n");

    /*
    pr->note[max_note_index].last_evidence_time = curr_time;
    */

    DEBUG(fprintf(stderr, "\t%g - hps %d\n", curr_time, pr->note[max_note_index].midi_number);)

    return;
}

/********************************************************************************************/
/* HSC_PITCH_ESTIMATE_2                                                                     */
/********************************************************************************************/
void hsc_pitch_estimate_2(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned char b, note, h, num_of_harm = 254, /* max_note_index = 0, */
    				num_of_notes = pr->num_of_notes;
    unsigned int harm_bin, next_harm_bin, num_of_bins = aud_mod->fft_p->num_of_bins;
    double hsc, band_hsc, /* max_val = (-1) * DBL_MAX, */ fund_pseudo_bin,
    		freq_res = aud_mod->fft_p->freq_resolution;

    for (note = 0; note < num_of_notes; note++)
    {
        if ( (fund_pseudo_bin = (pr->note[note].fund_freq / freq_res)) < 1 )
        	continue;

    	for (h = 1, hsc = band_hsc = b = 0; h <= num_of_harm; h++)
      	{
           harm_bin = (int)(fund_pseudo_bin * h - 0.5);
           next_harm_bin = (int)(fund_pseudo_bin * (h + 1) - 0.5);

           if (harm_bin >= num_of_bins)
           		break;

           if (h > pow(2, (double)b * 1/3)) /* otimizar com bitshift */
           {
               hsc += log10(band_hsc);
               band_hsc = 0;
               b++;
           }

           if (next_harm_bin < num_of_bins)
           		band_hsc += min(spec->content[harm_bin], spec->content[next_harm_bin]);
           else
           		band_hsc += spec->content[harm_bin];
        }
        if (band_hsc > 0)
        	hsc += log10(band_hsc);

     	printf("%g %u %g\n", curr_time, (note + 1), hsc);

     	/*
     	if (hsc > max_val)
     	{
     		max_val = hsc;
          	max_note_index = note;
    	}
    	*/
    }
    /* gnuplots needs this */
    printf("\n");

    /*
    pr->note[max_note_index].last_evidence_time = curr_time;
    */

    DEBUG(fprintf(stderr, "\t%g - hps %d\n", curr_time, pr->note[max_note_index].midi_number);)

    return;
}

/********************************************************************************************/
/* HSC_PITCH_ESTIMATE_3                                                                     */
/********************************************************************************************/
void hsc_pitch_estimate_3(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned char b, note, h, num_of_harm = 23, /* max_note_index = 0, */
    				num_of_notes = pr->num_of_notes,
        			first_note_midi_number = pr->note[0].midi_number;
    unsigned int harm_bin, next_harm_bin, num_of_bins = aud_mod->fft_p->num_of_bins;
    double hsc, band_hsc, /* max_val = (-1) * DBL_MAX, */ fund_pseudo_bin,
    freq_res = aud_mod->fft_p->freq_resolution, raw_freq_res = aud_mod->fft_p->freq_resolution
    													* aud_mod->fft_p->zero_padding_ratio;

    for (note = 0; note < num_of_notes; note++)
    {
        /*
    	printf("\n\n%g NOTE %u : %g hz ->", curr_time, note, pr->note[note].fund_freq);
    	*/

        if ( pr->note[note].fund_freq < raw_freq_res * 3.3 )
        	continue;

        fund_pseudo_bin = (pr->note[note].fund_freq / freq_res);

    	for (h = 1, hsc = band_hsc = b = 0; h <= num_of_harm; h++)
      	{
           harm_bin = (int)(fund_pseudo_bin * h - 0.5);
           next_harm_bin = (int)(fund_pseudo_bin * (h + 1) - 0.5);

           if (harm_bin >= num_of_bins)
           		break;

           		/*
           printf(" %u [%g]", harm_bin + 1, spec->content[harm_bin]);
           if (note == 50)
           		printf(" (24 [%g])", spec->content[23]);
           		*/

           if (h > pow(2, (double)b * 1/3)) /* otimizar com bitshift */
           {
               hsc += (band_hsc * (0.25 + 0.75 / (pow(2, b) - pow(2, b-1)) ) );
               band_hsc = 0;
               b++;
           }

           if (next_harm_bin < num_of_bins)
           		band_hsc += min(spec->content[harm_bin], spec->content[next_harm_bin]);
           else
           {
           		band_hsc += spec->content[harm_bin];
           		break;
           }
        }
        if (band_hsc > 0)
        	hsc += (band_hsc * (0.25 + 0.75 / (h-1 - pow(2, b-1)) ) );

        printf("%g %u %g\n", curr_time, (note + first_note_midi_number), hsc);

     	/*
     	if (hsc > max_val)
     	{
     		max_val = hsc;
          	max_note_index = note;
    	}
    	*/
    }
    /* gnuplots needs this */
    printf("\n");

    /*
    pr->note[max_note_index].last_evidence_time = curr_time;
    */

    DEBUG(fprintf(stderr, "\t%g - hps %d\n", curr_time, pr->note[max_note_index].midi_number);)

    return;
}


/********************************************************************************************/
/* BANDWISE_HSC_PITCH_ESTIMATE                                                              */
/********************************************************************************************/
/* expects POWER spectrum, maybe MAG, not LOG */
void bandwise_hsc_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    unsigned char note, b, h, num_of_harm = 7, most_salient_note_index = 0,
    				num_of_notes = pr->num_of_notes, num_of_bands = aud_mod->num_of_bands;
    unsigned int harm_bin, band_last_bin;
    double hsc, band_hsc, max_val = (-1) * DBL_MAX,
    		freq_res = aud_mod->fft_p->freq_resolution, fund_pseudo_bin, fund_true_bin;

    for (note = 0; note < num_of_notes; note++)
    {
        if (pr->note[note].fund_freq < freq_res)
        	continue;

    	/* ESTIMAÇÃO DE FATO */
    	fund_pseudo_bin = pr->note[note].fund_freq / freq_res;
    	fund_true_bin = fund_pseudo_bin - 1;
        hsc = 0;

        for (b = 0; b < num_of_bands; b++)
        {
            band_last_bin = aud_mod->band[b].first_bin + aud_mod->band[b].num_of_bins - 1;

            if (fund_true_bin > band_last_bin)
            	continue;

            band_hsc = 0;

        	for (h = (int)ceil((aud_mod->band[b].first_bin + 1) / fund_pseudo_bin);
         		 h < num_of_harm; h++)
          	{
                   harm_bin = (int)(fund_pseudo_bin * h - 0.5);

                   if (harm_bin > band_last_bin)
                   		break;

                   band_hsc += (spec->content[harm_bin]
                     				 * pow(10, aud_mod->band[b].response[harm_bin]));
            }

            if (band_hsc > 0)
            	hsc += log10(band_hsc);
        }

        printf("%g %u %g\n", curr_time, (note + 1), hsc);

     	if (hsc > max_val)
     	{
     		max_val = hsc;
          	most_salient_note_index = note;
    	}

    }
    /* gnuplots needs this */
    printf("\n");

    /*
    pr->note[max_note_index].last_evidence_time = curr_time;
    */

    DEBUG(fprintf(stderr, "\t%g - hps %d\n", curr_time, pr->note[max_note_index].midi_number);)

    return;
}


/********************************************************************************************/
/* FFT_FFT_PITCH_ESTIMATE                                                                   */
/********************************************************************************************/
void fft_fft_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time)
{
    static unsigned char first_call = TRUE;
    static Complex* v;

    if (spec == NULL || aud_mod == NULL || pr == NULL)
    {
        if (!first_call) free(v);
        return;
    }

    if (first_call)
    {
        v = (Complex*)emalloc(aud_mod->fft_p->num_of_bins * sizeof(Complex));

        first_call = FALSE;
    }
}

/********************************************************************************************/
/* CALC_SPECTRUM_OFFSET                                                                     */
/********************************************************************************************/
void calc_spectrum_offset(Buffer *buf, Spectrum_Type spec_type)
{
    void *p1 = &buf->spec_data[0][0], *p2;

    switch (spec_type)
    {
        case PH_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].phase);
            break;
        case UPH_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].unwrapped_phase);
            break;
        case MAG_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].magnitude);
            break;
        case POW_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].power);
            break;
        case LP_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].log_power);
            break;
        case WD_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].warped_denoised);
            break;
        case PN_SPECTRUM:
        	p2 = &(buf->spec_data[0][0].power_noise);
            break;
        default:
            fprintf(stderr, "Error: unknown spectrum type (code '%d').\n", spec_type);
            return;
    }
    CHOSEN_SPECTRUM_OFFSET = p2 - p1;
}

/********************************************************************************************/
/* CALC_SPECTRA_EXTRA_DATA                                                                  */
/********************************************************************************************/
void calc_spectra_extra_data(Spectra* spec, unsigned int num_of_bins)
{
    unsigned int k, max_index;
    double total_sum, max_val, curr_val;

    /* phase and noise spectra seem not to need the extra data, so they're not computed */

    for (k = 0, total_sum = max_index = 0, max_val = spec->magnitude.content[0];
    	 k < num_of_bins;
      	 k++)
   	{
   		total_sum += (curr_val = spec->magnitude.content[k]);
   		if (curr_val > max_val)
   		{
   			max_index = k;
   			max_val = curr_val;
     	}
   	}
   	spec->magnitude.max_index = max_index;
    spec->magnitude.total_sum = total_sum;

    for (k = 0, total_sum = max_index = 0, max_val = spec->power.content[0];
    	 k < num_of_bins;
      	 k++)
   	{
   		total_sum += (curr_val = spec->power.content[k]);
   		if (curr_val > max_val)
   		{
   			max_index = k;
   			max_val = curr_val;
     	}
   	}
   	spec->power.max_index = max_index;
    spec->power.total_sum = total_sum;

    for (k = 0, total_sum = max_index = 0, max_val = spec->log_power.content[0];
    	 k < num_of_bins;
      	 k++)
   	{
   		total_sum += (curr_val = spec->log_power.content[k]);
   		if (curr_val > max_val)
   		{
   			max_index = k;
   			max_val = curr_val;
     	}
   	}
   	spec->log_power.max_index = max_index;
    spec->log_power.total_sum = total_sum;

    for (k = 0, total_sum = max_index = 0, max_val = spec->warped_denoised.content[0];
    	 k < num_of_bins;
      	 k++)
   	{
   		total_sum += (curr_val = spec->warped_denoised.content[k]);
   		if (curr_val > max_val)
   		{
   			max_index = k;
   			max_val = curr_val;
     	}
   	}
   	spec->warped_denoised.max_index = max_index;
    spec->warped_denoised.total_sum = total_sum;
}


/********************************************************************************************/
/* SPL_TO_PHON                                                                              */
/********************************************************************************************/
double spl_to_phon(double pressure, double freq)
{
	const static unsigned char NUM_OF_SAMPLE_POINTS = 5, NUM_OF_CURVES = 10;
	static unsigned char first_call = TRUE;

	double curvesSamplePoints[] = {20, 100, 500, 2000, 4000, 9000, 12500, 16000};
	double **curvesValues;

	unsigned int i;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * I. Verify that frequency is within expected range                                     *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	if (freq < MIN_AUDIBLE_FREQUENCY  || freq > MAX_AUDIBLE_FREQUENCY) {
		return 0;
	} else if (freq < curvesSamplePoints[0]) {
		freq = curvesSamplePoints[0];
	} else if (freq > curvesSamplePoints[NUM_OF_SAMPLE_POINTS - 1]) {
		freq = curvesSamplePoints[NUM_OF_SAMPLE_POINTS - 1];
	}


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * II. Prepare for processing (initialize structures and parameters)                     *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	if (first_call) {
		curvesValues = (double**)emalloc(NUM_OF_SAMPLE_POINTS * sizeof(double*));
		for (i = 0; i < NUM_OF_SAMPLE_POINTS; i++) {
			curvesValues[i] = (double*)emalloc(NUM_OF_CURVES * sizeof(double*));
		}

		curvesValues[0][0] = -74;
		curvesValues[1][0] = -6;
		curvesValues[2][0] = 1;
		curvesValues[3][0] = 1;
		curvesValues[4][0] = 1;

		curvesValues[0][1] = 1;
		curvesValues[1][1] = 1;
		curvesValues[2][1] = 1;
		curvesValues[3][1] = 1;
		curvesValues[4][1] = 1;

		curvesValues[0][2] = 1;
		curvesValues[1][2] = 1;
		curvesValues[2][2] = 1;
		curvesValues[3][2] = 1;
		curvesValues[4][2] = 1;

		curvesValues[0][3] = 1;
		curvesValues[1][3] = 1;
		curvesValues[2][3] = 1;
		curvesValues[3][3] = 1;
		curvesValues[4][3] = 1;

		curvesValues[0][4] = 1;
		curvesValues[1][4] = 1;
		curvesValues[2][4] = 1;
		curvesValues[3][4] = 1;
		curvesValues[4][4] = 1;

		curvesValues[0][5] = 1;
		curvesValues[1][5] = 1;
		curvesValues[2][5] = 1;
		curvesValues[3][5] = 1;
		curvesValues[4][5] = 1;

		curvesValues[0][6] = 1;
		curvesValues[1][6] = 1;
		curvesValues[2][6] = 1;
		curvesValues[3][6] = 1;
		curvesValues[4][6] = 1;

		curvesValues[0][7] = 1;
		curvesValues[1][7] = 1;
		curvesValues[2][7] = 1;
		curvesValues[3][7] = 1;
		curvesValues[4][7] = 1;

		curvesValues[0][8] = 1;
		curvesValues[1][8] = 1;
		curvesValues[2][8] = 1;
		curvesValues[3][8] = 1;
		curvesValues[4][8] = 1;

		curvesValues[0][9] = 1;
		curvesValues[1][9] = 1;
		curvesValues[2][9] = 1;
		curvesValues[3][9] = 1;
		curvesValues[4][9] = 1;

		first_call = FALSE;
	}

	return 0;
}




















/********************************************************************************************/
/* FFT_PH_UPH_MAG_POW_LP                                                                    */
/********************************************************************************************/
void fft_ph_uph_mag_pow_lp(Buffer* buf, const unsigned char channel,
						const unsigned int sample_index, const unsigned int spec_index)
{
	static unsigned char first_call = TRUE;
	static double *windowed_signal = NULL;

	unsigned int i, core_length, total_length, last_bin;
	float *signal;
	double aux, cum, cum_incr, *fft_p_window_content,
			*spec_phase_content, *spec_unwrapped_phase_content, *spec_magnitude_content,
   			*spec_power_content, *spec_log_power_content;

	Spectra *spec;

	if (buf == NULL)
	{
        if (!first_call)
        {
        	if (windowed_signal) free(windowed_signal);
			finalize_fft_real_signal_to_power_phase();
        }
        return;
    }

    signal = &(buf->sample_data[channel][sample_index]);
	fft_p_window_content = &(buf->fft_p->window.content[0]);
	spec = &(buf->spec_data[channel][spec_index]);
	core_length = buf->fft_p->window.length;
	total_length = core_length * buf->fft_p->zero_padding_ratio;
	last_bin = buf->fft_p->num_of_bins - 1;

	if (first_call)
	{
    	windowed_signal = (double*) emalloc(total_length * sizeof(double));

		/* does the zero padding, if needed */
		for (i = core_length; i < total_length; i++)
			windowed_signal[i] = 0;

		prepare_fft_real_signal_to_power_phase(total_length);

		first_call = FALSE;
	}

   	for (i = 0; i < core_length; i++)
    		windowed_signal[i] = signal[i] * fft_p_window_content[i];

    spec_magnitude_content = spec->magnitude.content;
    spec_phase_content = spec->phase.content;
    spec_power_content = spec->power.content;
   	spec_log_power_content = spec->log_power.content;
   	spec_unwrapped_phase_content = spec->unwrapped_phase.content;

	fft_real_signal_to_power_phase(spec_power_content, spec_phase_content,
									windowed_signal, total_length);

  	for (i = 0, cum = cum_incr = (buf->first_sample_file_pos + sample_index)
  																	* TWO_PI / core_length;
   		 i <= last_bin;
       	 i++, cum += cum_incr)
    {
        spec_unwrapped_phase_content[i] = cum + spec_phase_content[i];

		aux = spec_power_content[i];
   		spec_power_content[i] = aux * POW_NORM_FACTOR;
   		spec_magnitude_content[i] = sqrt(aux) * MAG_NORM_FACTOR;
   		spec_log_power_content[i] = 10 * log10(aux) + LP_NORM_PARCEL;
    }

    calc_spectra_extra_data(spec, last_bin + 1);
}



/********************************************************************************************/
/* PRINT_SPEC_INTERP_TRIGONOMETRIC                                                          */
/********************************************************************************************/
void print_spec_interp_trigonometric(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double *spec_content,
  			t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;

  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].log_power.content;
  	for (k = 1; k < (num_of_bins - 1); k++)

      	if (spec_content[k - 1] <= spec_content[k] && spec_content[k + 1] <= spec_content[k])
      	{
      		const double accurate_bin = trigonometric_interp(buf, 0, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_bin - (double)k, &buf->fft_p->window);

      /*
			const double acc_left_bin = trigonometric_interp(buf, 0, spec_index, k-1),
						 acc_right_bin = trigonometric_interp(buf, 0, spec_index, k+1);

			double mean, var, max_var;

			mean = (acc_left_bin + accurate_bin + acc_right_bin)/3.0;
			var = ((acc_left_bin-mean)*(acc_left_bin-mean)
					+ (accurate_bin-mean)*(accurate_bin-mean)
					+ (acc_right_bin-mean)*(acc_right_bin-mean))/3.0;
			max_var = 1.0/3.0; */ /* also 1.0/6.0 */


      /* if (fabs(accurate_bin - k) <= 0.5 && var <= max_var) {
	 			printf("%g %g %g\n", t, freq_res * (1 + accurate_bin), accurate_log_power);
			}
      */
      printf("%g %g %g\n", t, freq_res * (1 + accurate_bin), accurate_log_power);

			/*
			if (fabs(accurate_bin - k) < 0.5) {
 				printf("%g %g %g\n", t, freq_res * (1 + accurate_bin), accurate_log_power);
 			}
			*/

 			/*
			if (fabs(accurate_bin - k) <= 0.52 && var < 1) {
     			printf("%g %g %g %g\n", t, freq_res * (1 + accurate_bin), accurate_log_power, var);
			}
			*/
      		/*
      		printf("%g %g %g %g\n", t, freq_res * (1+k), freq_res * (1 + accurate_freq), accurate_log_power);
      		*/
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
        spec_content = buf->spec_data[c][spec_index].log_power.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = trigonometric_interp(buf, c, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}


/********************************************************************************************/
/* PRINT_SPEC_INTERP_CHARPENTIER                                                            */
/********************************************************************************************/
void print_spec_interp_charpentier(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset)
{
  const double freq_res = buf->fft_p->freq_resolution;

  unsigned char c;
  unsigned int k, num_of_bins =  buf->fft_p->num_of_bins - 1;
  /* TODO: optimize with static variables or something similar to avoid this calculation */
  double *spec_content,
  			t = ((double)(buf->spec_to_sample_offset + (spec_index * buf->fft_p->step))
				+ buf->first_sample_file_pos) / buf->stream->sample_rate;


  /* TODO: possibly FILTER already here */
  if (buf->stream->num_of_channels == 1)
  {
  	/*
  	spec_content = ((Spectrum*)((void*)&(buf->spec_data[0][spec_index]) + offset))->content;
  	*/
  	spec_content = buf->spec_data[0][spec_index].log_power.content;
  	for (k = 1; k < (num_of_bins - 1); k++)

      	if (spec_content[k - 1] <= spec_content[k] && spec_content[k + 1] <= spec_content[k])
      	{
      		const double accurate_freq = charpentier_interp(buf, 0, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

      		printf("%g %g %g\n", t, freq_res * (1 + accurate_freq), accurate_log_power);
			k++; /* optimization based on the impossibility of two consecutives local maxima */
      	}
  }
  else
  {
      for (c = 0; c < buf->stream->num_of_channels; c++)
      {
      	/*
        spec_content = ((Spectrum*)((void*)&(buf->spec_data[c][spec_index]) + offset))->content;
        */
        spec_content = buf->spec_data[c][spec_index].log_power.content;
	  	for (k = 1; k < (num_of_bins - 1); k++)
      	if (spec_content[k - 1] < spec_content[k] && spec_content[k] > spec_content[k + 1])
   	  	{
      		const double accurate_freq = charpentier_interp(buf, c, spec_index, k),
				accurate_log_power = spec_content[k] - window_log_power_resp(accurate_freq - (double)k, &buf->fft_p->window);

        	printf("%g %u %g %g\n", t, c + 1, freq_res * (1 + accurate_freq), accurate_log_power);
        	k++; /* optimization based on the impossibility of two consecutives local maxima */
   	  	}
      }
  }
  printf("\n"); /* gnuplot needs this */
}

/********************************************************************************************/
/* TRIGONOMETRIC_INTERP                                                                     */
/********************************************************************************************/
double trigonometric_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp)
{
	const unsigned int core_len = buf->fft_p->window.length,
					total_len = core_len * buf->fft_p->zero_padding_ratio,
					num_of_bins = buf->fft_p->num_of_bins;
	const float *signal = &(buf->sample_data[channel][buf->fft_p->step * spec_index]);

	static unsigned char first_call = TRUE;
	static unsigned int last_spec_index = UINT_MAX;
	static signed long last_call_buf_first_sample_file_pos = LONG_MAX;
	static Complex *c_dh_win_spec = NULL;
	static double *windowed_signal = NULL;

	const double sample_rate = buf->stream->sample_rate,
			*power = buf->spec_data[channel][spec_index].power.content,
			*phase = buf->spec_data[channel][spec_index].phase.content;

	unsigned int i;
	double accurate_bin;

	if (buf == NULL)
    {
        if (!first_call)
        {
        	if (c_dh_win_spec) free(c_dh_win_spec);
        	if (windowed_signal) free(windowed_signal);

        	finalize_fft_real_signal();

        	first_call = TRUE;
        }

        return 0;
    }

	if (first_call == TRUE)
	{
		c_dh_win_spec = (Complex*)emalloc(total_len/2 * sizeof(Complex));
		windowed_signal = (double*)emalloc(total_len * sizeof(double));

		for (i = core_len; i < total_len; i++) {
			windowed_signal[i] = 0;
		}

		prepare_fft_real_signal(total_len);

		first_call = FALSE;
	}

	if (spec_index != last_spec_index
		|| last_call_buf_first_sample_file_pos != buf->first_sample_file_pos)
	{
		for (i = 0; i < core_len; i++) {
			windowed_signal[i] = (*(signal + i - 1) * buf->fft_p->window.content[i]);
		}
		fft_real_signal(c_dh_win_spec, windowed_signal, total_len);

		last_spec_index = spec_index;
		last_call_buf_first_sample_file_pos = buf->first_sample_file_pos;

		/*
		double time = (buf->first_sample_file_pos + buf->fft_p->step * spec_index)/(double)sample_rate,
				freq_res = buf->fft_p->freq_resolution;
		if (buf->first_sample_file_pos + buf->fft_p->step * spec_index >= 0) {
			for (i = 0; i < num_of_bins; i++) {
				fprintf(stderr, "%g %g %g\n", time, freq_res*(i+1), 10*log10(c_modulus(c_dh_win_spec[i])));
			}
			fprintf(stderr, "\n");
		}
		*/
	}

	Complex prev, curr;

	c_copy(c_polar_to_complex(sqrt(power[bin_to_interp]), phase[bin_to_interp]), &curr);
	c_copy(c_mult_by_real(c_dh_win_spec[bin_to_interp], MAG_NORM_FACTOR), &prev);

	if (bin_to_interp < num_of_bins/2) {
		accurate_bin  = (sample_rate*INVERSE_OF_PI * asin(
				c_modulus(c_subtract(curr, prev)) / (2*c_modulus(curr))
				))/ buf->fft_p->freq_resolution - 1;
	} else {
		accurate_bin  = (sample_rate*INVERSE_OF_PI * acos(
				c_modulus(c_add(curr, prev)) / (2*c_modulus(curr))
				))/ buf->fft_p->freq_resolution - 1;
	}

	return accurate_bin;
	/*
	if (fabs(accurate_bin - bin_to_interp) < 0.52)
		return accurate_bin;

	return bin_to_interp;
	*/
}
