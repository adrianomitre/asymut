#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "config_constants.h"
#include "config_variables.h"
#include "numeric_constants.h"
#include "aux_fun.h"
#include "spec_fun.h"
#include "stdft.h"
#include "buffer.h"
#include "notes.h"

/********************************************************************************************/
/* INIT_BUFFER                                                                              */
/********************************************************************************************/
signed char init_buffer(Buffer **buffer, Audio_Stream* stream, FFT_properties *fft_p)
{
	const unsigned int fft_size = fft_p->window.length * fft_p->zero_padding_ratio;
	unsigned int i, sample_array_length;
	unsigned int num_of_mov_median_elements, min_nomme, max_nomme;
	unsigned char chan;
	Buffer *buf;

	if (*buffer == NULL) *buffer = (Buffer*) emalloc(sizeof(Buffer));

	buf = *buffer;

    fft_p->freq_resolution = (double)stream->sample_rate / fft_size;

    fft_p->time_resolution = (double)fft_p->step / stream->sample_rate;

    min_nomme = ceil(ONSET_THRESHOLD_MIN_WINDOW_LENGTH/(fft_p->time_resolution));
	max_nomme = floor(ONSET_THRESHOLD_MAX_WINDOW_LENGTH/(fft_p->time_resolution));
	
	num_of_mov_median_elements = (min_nomme + max_nomme)/2;
	
	if (num_of_mov_median_elements % 2 == 0) {
		if (num_of_mov_median_elements > min_nomme) {
			num_of_mov_median_elements--;
		}
		else if (num_of_mov_median_elements < max_nomme) {
			num_of_mov_median_elements++;
		}
	}
	
	num_of_mov_median_elements = 1; /* TEMPORÁRIO, só até o esquema de mediana sair do Asymut */

    buf->threshold_to_unpredictability_offset = num_of_mov_median_elements/2;

    buf->sample_array_length = pow(2, ceil(log2(fft_p->window.length
    							+ fft_p->step * (num_of_mov_median_elements + 1))));

	sample_array_length = buf->sample_array_length;


	/* 1. Ensures that the parameters are consistent and all conditions are met */

	if (ONSET_THRESHOLD_MAX_WINDOW_LENGTH <= ONSET_THRESHOLD_MIN_WINDOW_LENGTH) {
		fprintf(stderr, "\nError: min_onset_win must be strictly smaller than max_onset_win.\n");
		return ERROR;
	}

	if (num_of_mov_median_elements % 2 == 0)
	{
		fprintf(stderr, "\n"
						"Error: there is no positive odd number of STDFT window steps\n"
						"       which satisfies 'min_onset_win' and 'max_onset_win'.\n"
  						"\n"
  						"Advice: try to enhance time resolution, by shortening the STDFT window\n"
  						"        length and/or hop, or enlarging the 'onset_win' range.\n");
		return ERROR;
	}
	if (!is_power_of_2(sample_array_length))
	{
		fprintf(stderr, "\nError: sample_array_length must be a power of two.\n");
		return ERROR;
	}
	if (sample_array_length < (fft_p->window.length + fft_p->step * (num_of_mov_median_elements + 1)))
	{
		fprintf(stderr, "\n"
						"Error: sample_array_length must be at least FFT_SIZE + STEP * (N_MEDIAN + 1).\n");
		return ERROR;
	}
	if (stream->total_num_of_samples < (fft_p->window.length + fft_p->step * (num_of_mov_median_elements - 1)))
	{
	    fprintf(stderr, "\nError: audio stream is of insufficient length.\n");
		return ERROR;
	}

	/* 2. Determine variable values */

	buf->spec_array_length = (sample_array_length - fft_p->window.length)/fft_p->step + 1;
 	buf->unpredictability_array_length = buf->spec_array_length - 2;
 	buf->threshold_array_length = buf->unpredictability_array_length - num_of_mov_median_elements + 1;

	buf->spec_to_sample_offset = fft_p->window.length/2;

	buf->unpredictability_to_spec_offset = 2;
 	buf->unpredictability_to_sample_offset = 2*fft_p->step + fft_p->window.length/2;

 	buf->threshold_win_len = num_of_mov_median_elements;
 	buf->threshold_to_unpredictability_offset = num_of_mov_median_elements/2;
 	buf->threshold_to_spec_offset = num_of_mov_median_elements/2 + 2;
 	buf->threshold_to_sample_offset = (num_of_mov_median_elements/2 + 2)*fft_p->step + fft_p->window.length/2;

  	buf->fft_p = fft_p;
 	buf->stream = stream;

	{
	 signed int tzl_cand = fft_p->window.length/2 - ((buf->stream->total_num_of_samples) % fft_p->step);

	 buf->trailing_zeroes_left = max(0, tzl_cand);
	}

 	buf->actual_num_of_specs = 0;
    buf->actual_num_of_unpredicts = 0;
    buf->actual_num_of_thresholds = 0;

 	/* 3. Allocates the structure elements */

  	buf->sample_data = (float**) emalloc(buf->stream->num_of_channels * sizeof(float*));
 	for (chan = 0; chan < buf->stream->num_of_channels; chan++) 	{
		buf->sample_data[chan] = (float*) emalloc(buf->sample_array_length * sizeof(float));
 	}

	if (OPERATION_MODE >= PRINT_SPEC) {
		buf->spec_data = (Spectra**) emalloc(buf->stream->num_of_channels * sizeof(Spectra*));
	 	for (chan = 0; chan < buf->stream->num_of_channels; chan++)
		{
			buf->spec_data[chan] = (Spectra*) emalloc(buf->spec_array_length * sizeof(Spectra));
			
			for (i = 0; i < buf->spec_array_length; i++)
			{
	        	buf->spec_data[chan][i].phase.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].unwrapped_phase.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].magnitude.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].power.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].log_power.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].warped_denoised.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
	        	buf->spec_data[chan][i].power_noise.content
	        							= (double*) emalloc(fft_p->num_of_bins * sizeof(double));
			}
		}
	}
	
	
	if (OPERATION_MODE >= PRINT_UNPRED) {
		buf->unpredictability_data = (double*) emalloc(buf->unpredictability_array_length
	 																		* sizeof(double));
	}
	if (OPERATION_MODE >= PRINT_THRESHOLD) {
		buf->threshold_data = (double*) emalloc(buf->threshold_array_length * sizeof(double));
	}

	/* 4. Does the initial padding with zeros */

	for (chan = 0; chan < buf->stream->num_of_channels; chan++)
		for (i = 0; i < buf->threshold_to_sample_offset; i++)
			buf->sample_data[chan][i] = 0;
			
	buf->actual_num_of_samples = buf->threshold_to_sample_offset;
    buf->first_sample_file_pos = (-1) * (signed long)buf->actual_num_of_samples;
    
  return OK;
}

/********************************************************************************************/
/* FILL_BUFFER                                                                              */
/********************************************************************************************/
signed char fill_buffer(Buffer* buf, Audio_Stream* input_stream, Auditory_Model* aud_mod,
							double unpredict_method(Buffer*, unsigned char, unsigned int))
{
	
    if (is_buffer_full_of_samples(buf)) {
        recycle_buffer(buf);
    }

    if (read_file_to_buffer(input_stream, buf) != OK)
    	return ERROR;

    if (is_buffer_full_of_samples(buf) == FALSE && buf->trailing_zeroes_left != 0) {
		pad_with_trailing_zeros(buf);
    }

    if (OPERATION_MODE > PRINT_SAMPLES)
    	return update_buffer(buf, aud_mod, unpredict_method);

   	return OK;
}

/********************************************************************************************/
/* IS_BUFFER_FULL_OF_SAMPLES                                                                */
/********************************************************************************************/
signed char is_buffer_full_of_samples(Buffer *buf)
{
	return (buf->actual_num_of_samples == buf->sample_array_length);
}

/********************************************************************************************/
/* IS_BUFFER_FULL_OF_THRESHOLD                                                              */
/********************************************************************************************/
signed char is_buffer_full_of_thresholds(Buffer *buf)
{
	return (buf->actual_num_of_thresholds == buf->threshold_array_length);
}

/********************************************************************************************/
/* READ_FILE_TO_BUFFER                                                                      */
/********************************************************************************************/
signed char read_file_to_buffer(Audio_Stream* input_stream, Buffer* buf)
{
    const unsigned char num_of_channels = input_stream->num_of_channels,
    		bytes_per_sample = input_stream->byte_depth * input_stream->num_of_channels;
    const unsigned int max_samples_to_read = buf->sample_array_length
    															- buf->actual_num_of_samples;
    unsigned int actually_read, last_index, i, j;
    unsigned char c, *raw_data;
    FILE *fp = input_stream->file_pointer;
    unsigned int const actual_num_of_samples = buf->actual_num_of_samples;

    if (max_samples_to_read == 0 || feof(fp) || ferror(fp) )
    	return ERROR;

	raw_data = (char*)emalloc(sizeof(char)*(max_samples_to_read * bytes_per_sample));

    actually_read = fread(raw_data, bytes_per_sample, max_samples_to_read, fp);

    last_index = (buf->actual_num_of_samples - 1) + actually_read;

    if (input_stream->bit_depth == 8)
    {
        for (c = 0; c < num_of_channels; c++)
            for (i = actual_num_of_samples, j = c;
          		 i <= last_index;
            	 i++, j += bytes_per_sample)
            {
            	buf->sample_data[c][i] = c8_to_float(raw_data[j]);
            }
    }
    else if (input_stream->bit_depth == 16)
    {
    	for (c = 0; c < num_of_channels; c++)
            for (i = actual_num_of_samples, j = 2*c;
            	 i <= last_index;
              	 i++, j += bytes_per_sample)
			{
            	buf->sample_data[c][i] = c16_to_float(&raw_data[j]);
			}
    }

    buf->actual_num_of_samples += actually_read;

    return OK;
}

/********************************************************************************************/
/* RECYCLE_BUFFER                                                                           */
/********************************************************************************************/
void recycle_buffer(Buffer* buf)
{
    unsigned char c;
    unsigned int i, first_to_keep, kept_number;

    /* recycle samples */

    first_to_keep = buf->threshold_array_length * buf->fft_p->step;
    kept_number = buf->sample_array_length - first_to_keep;

    for (c = 0; c < buf->stream->num_of_channels; c++)
        memmove(buf->sample_data[c],
                buf->sample_data[c] + first_to_keep,
                kept_number * sizeof(float));

    buf->actual_num_of_samples = kept_number;
    buf->first_sample_file_pos += first_to_keep;

    if (OPERATION_MODE <= PRINT_SAMPLES) {
    	return;
    }
    
    /* recycle spec */
    
    kept_number = 2 * buf->threshold_to_unpredictability_offset + 2;
    first_to_keep = buf->spec_array_length - kept_number;
    for (c = 0; c < buf->stream->num_of_channels; c++)
    	for (i = 0; i < kept_number; i++)
    		copy_spectra(buf, c, first_to_keep + i, i);

    buf->actual_num_of_specs = kept_number;
    
    if (OPERATION_MODE < PRINT_UNPRED) {
    	return;
    }
    
    /* recycle unpredictability */
    
    kept_number = 2 * buf->threshold_to_unpredictability_offset;
    first_to_keep = buf->unpredictability_array_length - kept_number;
    memmove(buf->unpredictability_data,
            buf->unpredictability_data + first_to_keep,
            kept_number * sizeof(double));
    buf->actual_num_of_unpredicts = kept_number;


    if (OPERATION_MODE <= PRINT_UNPRED) {
    	return;
    }

    /* 'recycle' thresholds */
    
    /* buf->recycled_threshold = buf->threshold_data[buf->actual_num_of_thresholds - 1]; */
    buf->actual_num_of_thresholds = 0;
}

/********************************************************************************************/
/* PAD_WITH_TRAILING_ZEROS                                                                  */
/********************************************************************************************/
void pad_with_trailing_zeros(Buffer* buf)
{
    unsigned int i;
    unsigned char c;

    if (buf->trailing_zeroes_left == 0)
    	return;

    if (buf->trailing_zeroes_left > buf->sample_array_length - buf->actual_num_of_samples)
    {
        for (c = 0; c < buf->stream->num_of_channels; c++)
        	for (i = buf->actual_num_of_samples; i < buf->sample_array_length; i++)
        		buf->sample_data[c][i] = 0;

        buf->trailing_zeroes_left -= (buf->sample_array_length - buf->actual_num_of_samples);
        buf->actual_num_of_samples = buf->sample_array_length;
    }
    else
    {
    	const unsigned int last_index = buf->actual_num_of_samples
    														+ buf->trailing_zeroes_left - 1;
		
        for (c = 0; c < buf->stream->num_of_channels; c++)
        	for (i = buf->actual_num_of_samples; i <= last_index; i++)
        		buf->sample_data[c][i] = 0;

        buf->trailing_zeroes_left = 0;
        buf->actual_num_of_samples = last_index + 1;
    }
}

/********************************************************************************************/
/* UPDATE_BUFFER                                                                            */
/********************************************************************************************/
signed char update_buffer(Buffer* buf, Auditory_Model* aud_mod,
 								double unpredict_method(Buffer*, unsigned char, unsigned int))
{
    const unsigned char num_of_channels = buf->stream->num_of_channels;
    const unsigned int fft_p_step = buf->fft_p->step;
    
    unsigned char c;
    unsigned int i, last_index;
    unsigned int actual_num_of_specs = buf->actual_num_of_specs;
    double *unpredictability_data = buf->unpredictability_data;

	if (
		is_buffer_full_of_thresholds(buf)
		 ||
		buf->actual_num_of_samples
		 <
		(buf->threshold_to_sample_offset + fft_p_step*(buf->threshold_to_unpredictability_offset + 1))

		/*
		buf->actual_num_of_samples < buf->threshold_to_sample_offset + buf->fft_p->window.length/2
		buf->actual_num_of_samples < 2*(buf->threshold_to_sample_offset - fft_p_step)
		*/
	   )
	{
		return ERROR;
	}

/*
    if (buf->actual_num_of_samples >= 2*(buf->threshold_to_sample_offset - fft_p_step)
        && is_buffer_full_of_thresholds(buf) == FALSE)
    {
 */ 
 
	/* I. Updates the Spectra (FFT) layer */
	
	last_index = (buf->actual_num_of_samples - buf->fft_p->window.length)
																/ fft_p_step;
	for (c = 0; c < num_of_channels; c++)
	    for (i = actual_num_of_specs; i <= last_index; i++)
	 	{
	    	fft_ph_uph_mag_pow_lp(buf, c, i * fft_p_step, i);
		}
	buf->actual_num_of_specs = actual_num_of_specs = last_index + 1;

	/* II. Updates the Unpredictability layer */
	
	if (OPERATION_MODE >= PRINT_UNPRED)
	{
	    last_index = actual_num_of_specs - 3;
	    for (i = buf->actual_num_of_unpredicts; i <= last_index; i++)
	    	for (c = 0, unpredictability_data[i] = 0; c < num_of_channels; c++)
	    	{
	         	/* TODO: what about a sophisticated criteria for combining
	          	unpredictability across channels? */
	
	    		unpredictability_data[i] += unpredict_method(buf, c, i);
	      										/*
	                							incr_log_power_unpredict(buf, c, i)
	                								* pow(phase_unpredict(buf, c, i), 16);
	
	                							*/
	    	}
	    buf->actual_num_of_unpredicts = last_index + 1;
	}
	
	if (OPERATION_MODE >= PRINT_THRESHOLD)
	{
	    calc_threshold(buf, (unsigned short)(buf->threshold_win_len
	    	* ONSET_THRESHOLD_PERCENTILE - 0.5), ONSET_THRESHOLD_SCALING_FACTOR,
	     	ONSET_THRESHOLD_CONST_PART);
	
	    buf->actual_num_of_thresholds = buf->actual_num_of_unpredicts
	 													- buf->threshold_win_len + 1;
	}		

	return OK;
}

/********************************************************************************************/
/* CALC_THRESHOLD                                                                           */
/********************************************************************************************/
void calc_threshold(Buffer* buf, unsigned short med_index, double scale_factor,
																			double const_part)
{
 	calc_mov_median(buf->unpredictability_data, buf->actual_num_of_unpredicts,
    				buf->threshold_win_len, med_index, scale_factor, const_part,
        			buf->threshold_data);
}


/********************************************************************************************/
/* TRANSCRIBE_BUF                                                                           */
/********************************************************************************************/
void transcribe_buf(Buffer* buf, Auditory_Model* aud_mod, Pitch_Range* pr,
					void (*pitch_estimate)(Spectrum*, Auditory_Model*, Pitch_Range*, double))
{
	static unsigned char flag = DOWN;
	static double last_onset_time = (-1) * DBL_MAX, last_flag_down_time = 0;
	
 	unsigned char ch;
 	unsigned int spec_index, unp_index, thr_index;
	double curr_time, time_step = (double) buf->fft_p->step / buf->stream->sample_rate;

	for (thr_index = 0,
 		 spec_index = buf->threshold_to_spec_offset,
     	 unp_index = buf->threshold_to_unpredictability_offset,
     	 curr_time = (double)(buf->threshold_to_sample_offset + buf->first_sample_file_pos)
				     	 / buf->stream->sample_rate;

    	 thr_index < buf->actual_num_of_thresholds;

    	 thr_index++,
      	 spec_index++,
         unp_index++,
         curr_time += time_step)
	{
	   	if (buf->unpredictability_data[unp_index] > buf->threshold_data[thr_index])
	   	{
     		if (flag == DOWN)
    	   	{
             	last_onset_time = curr_time;
             	flag = UP;
            }
        }
        else if ((curr_time - last_flag_down_time) >= MIN_INTERONSET_DISTANCE)
        {
        	flag = DOWN;
        	last_flag_down_time = curr_time;
        }

        if (pr->current_polyphony == 0
        	&& (curr_time - last_onset_time) > MAX_DELAY_AFTER_ONSET)
        	continue;
        	
        for (ch = 0; ch < buf->stream->num_of_channels; ch++)
        {
       		pitch_estimate((Spectrum*)((void*)&buf->spec_data[ch][spec_index]
         					+ CHOSEN_SPECTRUM_OFFSET), aud_mod, pr, curr_time);
        }
        
        if ((curr_time - last_onset_time) <= MAX_DELAY_AFTER_ONSET)
			turn_on_new_notes(pr, last_onset_time, curr_time);
			
        turn_off_obsolete_notes(pr, curr_time);
	}

	if (buf->trailing_zeroes_left == 0)
	{
        turn_off_obsolete_notes(pr, DBL_MAX);
        if (PRINT_NOTE_AS_GNUPLOT_ARROWS)
        {
        	printf("set yrange [%u:%u]\n", pr->note[0].midi_number,
		        								pr->note[0].midi_number + pr->num_of_notes);
		    printf("set xrange [0:%g]\n", curr_time);
		    printf("plot 0\n");
        }
	}
}

/********************************************************************************************/
/* FREE_BUFFER                                                                              */
/********************************************************************************************/
void free_buffer(Buffer* buf)
{
    unsigned char c;
    unsigned int i;

    if (buf->sample_array_length == 0)
    	return;

    for (c = 0; c < buf->stream->num_of_channels; c++)
    {
    	free(buf->sample_data[c]);
    	if (OPERATION_MODE >= PRINT_SPEC) {
	    	for (i = 0; i < buf->spec_array_length; i++)
			{
	        	free(buf->spec_data[c][i].phase.content);
	        	free(buf->spec_data[c][i].unwrapped_phase.content);
	        	free(buf->spec_data[c][i].magnitude.content);
	        	free(buf->spec_data[c][i].power.content);
	        	free(buf->spec_data[c][i].log_power.content);
	        	free(buf->spec_data[c][i].warped_denoised.content);
	        	free(buf->spec_data[c][i].power_noise.content);
			}
	    	free(buf->spec_data[c]);
    	}
    }
    free(buf->unpredictability_data);
    free(buf->threshold_data);

    buf->sample_array_length = 0;
    
    free_window(&buf->fft_p->window);

    fft_ph_uph_mag_pow_lp(NULL, 0, 0, 0);
}

/********************************************************************************************/
/* IS_BUF_CONSISTENT                                                                        */
/********************************************************************************************/
signed char is_buffer_consistent(Buffer* buf, FFT_properties* fft_p)
{
	if (buf->actual_num_of_samples > buf->sample_array_length)
	{
		fprintf(stderr, "Buffer inconsistency: 'actual_num_of_samples' value.\n");
		return ERROR;
    }
	if (buf->actual_num_of_specs > ((buf->actual_num_of_samples - fft_p->window.length)
										/fft_p->step + 1))
	{
		fprintf(stderr, "Buffer inconsistency: actual_num_of_specs' value.\n");
		return ERROR;
	}
	if (buf->actual_num_of_unpredicts > buf->actual_num_of_specs - 2)
	{
		fprintf(stderr, "Buffer inconsistency: 'actual_num_of_unpredicts' value.\n");
		return ERROR;
	}
	if (buf->actual_num_of_thresholds > (buf->actual_num_of_unpredicts
 											+ buf->threshold_array_length
            								- buf->unpredictability_array_length))
	{
		fprintf(stderr, "Buffer inconsistency: 'actual_num_of_thresholds' value.\n");
		return ERROR;
	}
	return OK;
}

/********************************************************************************************/
/* PRINT_BUFFER_STATE                                                                       */
/********************************************************************************************/
void print_buffer_state(Buffer* buf, const int n)
{
    fprintf(stdout, "***********************\n");
    fprintf(stdout, "****** PBS %d *********\n", n);
    fprintf(stdout, "***********************\n");
    fprintf(stdout, "threshold_win_len: %u\n", buf->threshold_win_len);
    fprintf(stdout, "buf->stream->num_of_channels: %u\n",  buf->stream->num_of_channels);
    fprintf(stdout, "sample_array_length: %u\n",  buf->sample_array_length);
    fprintf(stdout, "spec_to_sample_offset %u\n", buf->spec_to_sample_offset);
    fprintf(stdout, "spec_array_length %u\n", buf->spec_array_length);
    fprintf(stdout, "unpredictability_to_sample_offset %u\n",
    												buf->unpredictability_to_sample_offset);
    fprintf(stdout, "unpredictability_to_spec_offset %u\n",
    												buf->unpredictability_to_spec_offset);
    fprintf(stdout, "unpredictability_array_length %u\n", buf->unpredictability_array_length);
    fprintf(stdout, "threshold_to_sample_offset %u\n", buf->threshold_to_sample_offset);
    fprintf(stdout, "threshold_to_spec_offset %u\n", buf->threshold_to_spec_offset);
    fprintf(stdout, "threshold_to_unpredictability_offset %u\n",
    											buf->threshold_to_unpredictability_offset);
    fprintf(stdout, "threshold_array_length %u\n", buf->threshold_array_length);
    fprintf(stdout, "actual_num_of_samples: %u\n",  buf->actual_num_of_samples);
    fprintf(stdout, "actual_num_of_specs %u\n", buf->actual_num_of_specs);
    fprintf(stdout, "actual_num_of_unpredicts %u\n", buf->actual_num_of_unpredicts);
    fprintf(stdout, "actual_num_of_thresholds %u\n", buf->actual_num_of_thresholds);
    fprintf(stdout, "trailing_zeroes_left %u\n", buf->trailing_zeroes_left);
    fprintf(stdout, "first_sample_file_pos %ld\n", buf->first_sample_file_pos);
}

/********************************************************************************************/
/* PRINT_BUF_SAMPLES                                                                        */
/********************************************************************************************/
void print_buf_samples(Buffer* buf)
{
    unsigned char c;
    unsigned int i;
	double t = buf->first_sample_file_pos / (double)buf->stream->sample_rate,
 			dt = 1 / (double)buf->stream->sample_rate;

 	unsigned int first_index, last_index = buf->threshold_array_length * buf->fft_p->step - 1;

 	if (buf->first_sample_file_pos < 0)
 		first_index = (-1) * buf->first_sample_file_pos;
 	else
 		first_index = 0;

	if (buf->stream->num_of_channels == 1)
	{
     	for (i = first_index; i <= last_index; i++)
     		printf("%g %g\n", t + i * dt, buf->sample_data[0][i]);
	}
	else
	{
        for (c = 0; c < buf->stream->num_of_channels; c++, t += dt)
         	for (i = first_index; i <= last_index; i++)
        		printf("%g %u %g\n", t + i * dt, c + 1, buf->sample_data[c][i]);

	}

    if (buf->trailing_zeroes_left == 0)
	{
    	if (buf->stream->num_of_channels == 1)
    	{
         	for (i = last_index + 1; i < buf->actual_num_of_samples; i++)
         		printf("%g %g\n", t + i * dt, buf->sample_data[0][i]);
    	}
    	else
    	{
            for (c = 0; c < buf->stream->num_of_channels; c++)
            	for (i = last_index + 1; i < buf->actual_num_of_samples; i++)
            		printf("%g %u %g\n", t + i * dt, c + 1, buf->sample_data[c][i]);

    	}
    }
}

/********************************************************************************************/
/* PRINT_BUF_INDEP_F0_ESTIMATES                                                             */
/********************************************************************************************/
void print_buf_indep_f0_estimates(Buffer* buf,
					void (*f0_method)(const float*, const unsigned int, double*, double*))
{
	const double dt = buf->fft_p->step / (double)buf->stream->sample_rate;
	unsigned int last_index;

 	unsigned int i;
 	int offset;
	double t;

	offset = buf->unpredictability_to_sample_offset - buf->fft_p->window.length/2;
	t = (offset + buf->first_sample_file_pos + buf->fft_p->window.length/2) / (double)buf->stream->sample_rate;
	last_index =  floor((buf->actual_num_of_samples - offset - buf->fft_p->window.length)/(double)buf->fft_p->step);
	
	/*
	printf("step: %d\n", buf->fft_p->step);
	printf("last_index: %d\n", last_index);
	*/

	if (buf->stream->num_of_channels == 1)
	{
		double period_estimate, estimate_quality;
		
     	for (i = 0; i <= last_index; i++) {
	     	f0_method(buf->sample_data[0] + offset + i*buf->fft_p->step, buf->fft_p->window.length,
	     					&period_estimate, &estimate_quality);
			printf("%g %g %g\n", t + i*dt, buf->stream->sample_rate/period_estimate, estimate_quality);
     	}
	}
	else
	{
		printf("ERROR: print_buf_indep_f0_estimates() does is not yet tuned for stereo files.");
		/*
        for (c = 0; c < buf->stream->num_of_channels; c++, t += dt)
         	for (i = 0; i <= last_index; i++)
        		printf("%g %u %g\n", t + i * dt, c + 1, buf->sample_data[c][i]);
		*/
	}
}

/********************************************************************************************/
/* PRINT_BUF_SPEC                                                                           */
/********************************************************************************************/
void print_buf_spec(Buffer* buf, void (*print_spec)(Buffer*, const unsigned int, const unsigned int))
{
	unsigned int i;
	unsigned int specs_to_print = buf->actual_num_of_specs;

	if (specs_to_print == buf->spec_array_length)
		specs_to_print -= buf->threshold_to_unpredictability_offset;

    for (i = buf->threshold_to_spec_offset; i < specs_to_print; i++)
  		print_spec(buf, i, CHOSEN_SPECTRUM_OFFSET);
}

/********************************************************************************************/
/* PRINT_BUF_UNPRED                                                                         */
/********************************************************************************************/
void print_buf_unpred(Buffer* buf)
{
	unsigned int i = buf->threshold_to_unpredictability_offset;
	const unsigned int last_index = min(buf->actual_num_of_unpredicts,
    				 					buf->unpredictability_array_length - i - 1);

    for (; i <= last_index; i++)
  		print_unpred(buf, i);
}

/********************************************************************************************/
/* PRINT_BUF_THRESHOLD                                                                      */
/********************************************************************************************/
void print_buf_threshold(Buffer* buf)
{
	unsigned int i, actual_num_of_thresholds = buf->actual_num_of_thresholds;

    for (i = 0; i < actual_num_of_thresholds; i++)
    	print_threshold(buf, i);
}
	
/********************************************************************************************/
/* PRINT_BUF_PITCH_ESTIMATE                                                                 */
/********************************************************************************************/
void print_buf_pitch_estimate(Buffer* buf, Auditory_Model* aud_mod, Pitch_Range* pr,
					void (*pitch_estimate)(Spectrum*, Auditory_Model*, Pitch_Range*, double))
{
 	unsigned char ch;
  	unsigned char const num_of_channels = buf->stream->num_of_channels;
 	unsigned int spec_index, unp_index, thr_index;
	double curr_time, time_step = (double) buf->fft_p->step / buf->stream->sample_rate;

	for (thr_index = 0,
 		 spec_index = buf->threshold_to_spec_offset,
     	 unp_index = buf->threshold_to_unpredictability_offset,
     	 curr_time = (double)(buf->threshold_to_sample_offset + buf->first_sample_file_pos)
				     	 / buf->stream->sample_rate;

    	 thr_index < buf->actual_num_of_thresholds;

    	 thr_index++,
      	 spec_index++,
         unp_index++,
         curr_time += time_step)
	{
        for (ch = 0; ch < num_of_channels; ch++)
       		pitch_estimate((Spectrum*)((void*)&buf->spec_data[ch][spec_index]
         					+ CHOSEN_SPECTRUM_OFFSET), aud_mod, pr, curr_time);
	}
}

/********************************************************************************************/
/* PRINT_BUF_ONSET                                                                          */
/********************************************************************************************/
void print_buf_onset(Buffer* buf)
{
	static unsigned char flag = DOWN;
	static double last_onset_time = (-1) * DBL_MAX, last_flag_down_time = 0;
	
	unsigned int spec_index, unp_index, thr_index;
	const unsigned int actual_num_of_thresholds = buf->actual_num_of_thresholds;

	double curr_time;
	const double *unpredictability_data = buf->unpredictability_data,
					*threshold_data = buf->threshold_data;
	const double time_step = (double) buf->fft_p->step / buf->stream->sample_rate;
	

	for (thr_index = 0,
 		 spec_index = buf->threshold_to_spec_offset,
     	 unp_index = buf->threshold_to_unpredictability_offset,
     	 curr_time = (double)(buf->threshold_to_sample_offset + buf->first_sample_file_pos)
				     	 / buf->stream->sample_rate;

    	 thr_index < actual_num_of_thresholds;

    	 thr_index++,
      	 spec_index++,
         unp_index++,
         curr_time += time_step)
	{
		if (flag == DOWN)
		{
    	   	if (unpredictability_data[unp_index] > threshold_data[thr_index]
         		&& (curr_time - last_flag_down_time) >= MIN_INTERONSET_DISTANCE)
        	   	{
                 	last_onset_time = curr_time;
                 	flag = UP;
                 	printf("%g\n", curr_time);
                }
        }
        else
        {
            if (unpredictability_data[unp_index] < threshold_data[thr_index])
            {
            	flag = DOWN;
            	last_flag_down_time = curr_time;
            }
        }
	}
}
