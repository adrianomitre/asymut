#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "numeric_constants.h"
#include "trigonometric_constants.h"
#include "config_constants.h"
#include "config_variables.h"
#include "stdft.h"
#include "complex.h"
#include "aux_fun.h"
#include "spec_fun.h"
#include "onset_detection.h"


/********************************************************************************************/
/* PARSE_ARGUMENTS                                                                          */
/********************************************************************************************/
signed char parse_arguments(const int n_args, char** argument,
                        	Audio_Stream* input, FFT_properties* fft)
{
    signed char status = OK;
    unsigned char enough_args = 0, step_set = FALSE, spec_type_set = FALSE;
    char wav_ext[5];
    int i;

    strcpy(wav_ext, ".wav"); /* TODO: verify if THIS is indeed necessary or omittable */

	/* set default values */

	PITCH_RANGE_LOWEST_NOTE = DEFAULT_PITCH_RANGE_LOWEST_NOTE;
	PITCH_RANGE_HIGHEST_NOTE = DEFAULT_PITCH_RANGE_HIGHEST_NOTE;
    MIN_ABSOLUTE_F0 = DEFAULT_MIN_ABSOLUTE_F0;
    MIN_F0_TO_FREQ_RES_RATIO = DEFAULT_MIN_F0_TO_FREQ_RES_RATIO;
    MAX_ABSOLUTE_F0 = DEFAULT_MAX_ABSOLUTE_F0;
    MAX_CRITICAL_BANDS = DEFAULT_MAX_CRITICAL_BANDS;
    MIN_INTERONSET_DISTANCE = DEFAULT_MIN_INTERONSET_DISTANCE;
    ONSET_THRESHOLD_MIN_WINDOW_LENGTH = DEFAULT_ONSET_THRESHOLD_MIN_WINDOW_LENGTH;
    ONSET_THRESHOLD_MAX_WINDOW_LENGTH = DEFAULT_ONSET_THRESHOLD_MAX_WINDOW_LENGTH;
    ONSET_THRESHOLD_PERCENTILE = DEFAULT_ONSET_THRESHOLD_PERCENTILE;
    ONSET_THRESHOLD_SCALING_FACTOR = DEFAULT_ONSET_THRESHOLD_SCALING_FACTOR;
    ONSET_THRESHOLD_CONST_PART = DEFAULT_ONSET_THRESHOLD_CONST_PART;
    MAX_DELAY_AFTER_ONSET = DEFAULT_MAX_DELAY_AFTER_ONSET;
    MAX_GAP_INSIDE_NOTE = DEFAULT_MAX_GAP_INSIDE_NOTE;
    MIN_NOTE_DURATION = DEFAULT_MIN_NOTE_DURATION;
    FREQ_REF_A4 = DEFAULT_FREQ_REF_A4;
    KLAPURI_ITERATION_CONTROL = DEFAULT_KLAPURI_ITERATION_CONTROL;
	fft->window.apodization_function.type = DEFAULT_APODIZATION_FUNCTION_TYPE;
    ESTIMATION_METHOD = DEFAULT_ESTIMATION_METHOD;
    UNPRED_METHOD = DEFAULT_UNPRED_METHOD;
    OPERATION_MODE = DEFAULT_OPERATION_MODE;
    NORMALIZATION_LEVEL = DEFAULT_NORMALIZATION_LEVEL;
    NORMALIZATION_ENABLED = DEFAULT_NORMALIZATION_ENABLED;
    DOUBLED_LINEBREAKS = FALSE;
    PRINT_NOTE_AS_GNUPLOT_ARROWS = FALSE;
    fft->zero_padding_ratio = DEFAULT_ZERO_PADDING_RATIO;
    fft->step = 0;

	/* parse arguments */

	for (i = 1; i < n_args; i++)
	{
		if (!strcmp(argument[i], "--fft_win_size") || !strcmp(argument[i], "-w"))
		{
      		if (++i < n_args)
        	{
             	if (strchr(argument[i], 'k')) fft->window.length = atof(argument[i]) * 1024;
             	else fft->window.length = atoi(argument[i]);
                enough_args++;
        	}
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--fft_win_step") || !strcmp(argument[i], "-s"))
		{
      		if (++i < n_args)
        	{
             	if (strchr(argument[i], 'k')) fft->step = atof(argument[i]) * 1024;
             	else fft->step = atoi(argument[i]);
             	step_set = TRUE;
            }
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--const_part") || !strcmp(argument[i], "-c"))
		{
      		if (++i < n_args)
             	ONSET_THRESHOLD_CONST_PART = atof(argument[i]);
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--zero_padding") || !strcmp(argument[i], "-z"))
		{
      		if (++i < n_args) fft->zero_padding_ratio = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--input_file") || !strcmp(argument[i], "-i"))
		{
      		if (++i < n_args)
        	{
             	strncpy(input->file_name, argument[i], FILENAME_MAX);
             	enough_args++;
        	}
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--min_abs_f0"))
		{
      		if (++i < n_args) MIN_ABSOLUTE_F0 = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--min_rel_f0"))
		{
      		if (++i < n_args) MIN_F0_TO_FREQ_RES_RATIO = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--max_abs_f0"))
		{
      		if (++i < n_args) MAX_ABSOLUTE_F0 = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--max_bands"))
		{
      		if (++i < n_args) MAX_CRITICAL_BANDS = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--min_interonset_dist"))
		{
      		if (++i < n_args) MIN_INTERONSET_DISTANCE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--min_onset_win"))
		{
      		if (++i < n_args) ONSET_THRESHOLD_MIN_WINDOW_LENGTH = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--max_onset_win"))
		{
      		if (++i < n_args) ONSET_THRESHOLD_MAX_WINDOW_LENGTH = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--percentile") || !strcmp(argument[i], "-p"))
		{
      		if (++i < n_args) ONSET_THRESHOLD_PERCENTILE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--scaling_factor") || !strcmp(argument[i], "-f"))
		{
      		if (++i < n_args) ONSET_THRESHOLD_SCALING_FACTOR = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--max_delay"))
		{
      		if (++i < n_args) MAX_DELAY_AFTER_ONSET = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--max_gap"))
		{
      		if (++i < n_args) MAX_GAP_INSIDE_NOTE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--min_duration"))
		{
      		if (++i < n_args) MIN_NOTE_DURATION = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--ref_a4_freq"))
		{
      		if (++i < n_args) FREQ_REF_A4 = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--pow_norm_level")  || !strcmp(argument[i], "-n"))
		{
      		if (++i < n_args) NORMALIZATION_LEVEL = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--lowest_note"))
        {
      		if (++i < n_args) PITCH_RANGE_LOWEST_NOTE = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--highest_note"))
        {
      		if (++i < n_args) PITCH_RANGE_HIGHEST_NOTE = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--f0_method"))
		{
      		if (++i < n_args)
      		{
            	if (!strcmp(argument[i], "hps")) ESTIMATION_METHOD = HPS;
             	else if (!strcmp(argument[i], "hsc")) ESTIMATION_METHOD = HSC;
             	else if (!strcmp(argument[i], "bw_hsc")) ESTIMATION_METHOD = BANDWISE_HSC;
             	else if (!strcmp(argument[i], "max_index")) ESTIMATION_METHOD = MAX_INDEX;
                else if (!strcmp(argument[i], "fft_fft")) ESTIMATION_METHOD = FFT_FFT;
                else if (!strcmp(argument[i], "klapuri")) ESTIMATION_METHOD = KLAPURI;
                else
                {
                    fprintf(stderr, "Unrecognized F0 estimation method: '%s'\n", argument[i]);
                    status = ERROR;
                }
                if (status != ERROR && ++i < n_args)
        	    {
                    spec_type_set = TRUE;
                    /* it doesn't make sense to use phase or noise spectrum here */
                    if (!strcmp(argument[i], "mag")) CHOSEN_SPECTRUM = MAG_SPECTRUM;
                    else if (!strcmp(argument[i], "pow")) CHOSEN_SPECTRUM = POW_SPECTRUM;
                    else if (!strcmp(argument[i], "lp")) CHOSEN_SPECTRUM = LP_SPECTRUM;
                    else if (!strcmp(argument[i], "wd")) CHOSEN_SPECTRUM = WD_SPECTRUM;
                    else if (!strcmp(argument[i], "ph")
                    		|| !strcmp(argument[i], "uph")
                    		|| !strcmp(argument[i], "pn"))
                    {
                    	fprintf(stderr, "Phase or noise spectrum are not useful in F0 "
                     					"estimation. Choose another type.\n");
                        spec_type_set = FALSE;
                    	status = ERROR;
                    }
                    else
                    {
                        spec_type_set = FALSE;
                        i--;
                    }
                }
                else
                {
                    switch (ESTIMATION_METHOD)
                    {
                        case MAX_INDEX:
                         	CHOSEN_SPECTRUM = DEFAULT_MAX_INDEX_SPECTRUM;
                            break;
                        case HPS:
                         	CHOSEN_SPECTRUM = DEFAULT_HPS_SPECTRUM;
                            break;
                        case HSC:
                         	CHOSEN_SPECTRUM = DEFAULT_HSC_SPECTRUM;
                            break;
                        case BANDWISE_HSC:
                         	CHOSEN_SPECTRUM = DEFAULT_BANDWISE_HSC_SPECTRUM;
                            break;
                        case FFT_FFT:
                         	CHOSEN_SPECTRUM = DEFAULT_FFT_FFT_SPECTRUM;
                            break;
                        case KLAPURI:
                         	CHOSEN_SPECTRUM = DEFAULT_KLAPURI_SPECTRUM;
                            break;
                    }
                }
            }
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--unpred_method"))
		{
      		if (++i < n_args)
      		{
            	if (!strcmp(argument[i], "lp")) UNPRED_METHOD = LP_UM;
             	else if (!strcmp(argument[i], "ilp")) UNPRED_METHOD = ILP_UM;
             	else if (!strcmp(argument[i], "ph")) UNPRED_METHOD = PH_UM;
             	else if (!strcmp(argument[i], "mph")) UNPRED_METHOD = MPH_UM;
                else if (!strcmp(argument[i], "complex")) UNPRED_METHOD = COMPLEX_UM;
                else if (!strcmp(argument[i], "nd_complex")) UNPRED_METHOD = NON_DECREASING_COMPLEX_UM;
                else if (!strcmp(argument[i], "new")) UNPRED_METHOD = NEW_UM;
                else if (!strcmp(argument[i], "rms")) UNPRED_METHOD = RMS_UM;
                else if (!strcmp(argument[i], "wrms")) UNPRED_METHOD = WRMS_UM;
                else if (!strcmp(argument[i], "erb")) UNPRED_METHOD = ERB_UM;
                else if (!strcmp(argument[i], "bark")) UNPRED_METHOD = BARK_UM;
                else if (!strcmp(argument[i], "mkl")) UNPRED_METHOD = MKL_UM;
             	else if (!strcmp(argument[i], "sd")) UNPRED_METHOD = SPEC_DIF_UM;
                else
                {
                    fprintf(stderr, "Unrecognized F0 estimation method: '%s'\n", argument[i]);
                    status = ERROR;
                }
            }
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--fft_apod_fun")  || !strcmp(argument[i], "-a"))
		{
      		if (++i < n_args)
      		{
            	if (!strcmp(argument[i], "rectangular")) fft->window.apodization_function.type = RECTANGULAR;
             	else if (!strcmp(argument[i], "triangular")) fft->window.apodization_function.type = TRIANGULAR;
             	else if (!strcmp(argument[i], "hamming")) fft->window.apodization_function.type = HAMMING;
                else if (!strcmp(argument[i], "hann")) fft->window.apodization_function.type = HANN;
                else if (!strcmp(argument[i], "blackman")) fft->window.apodization_function.type = BLACKMAN;
                else if (!strcmp(argument[i], "blackman-harris")) fft->window.apodization_function.type = BLACKMAN_HARRIS;
                else if (!strcmp(argument[i], "nuttall3")) fft->window.apodization_function.type = NUTTALL3;
                else if (!strcmp(argument[i], "nuttall4")) fft->window.apodization_function.type = NUTTALL4;
                else if (!strcmp(argument[i], "nuttall5")) fft->window.apodization_function.type = NUTTALL5;
                else if (!strcmp(argument[i], "nuttall6")) fft->window.apodization_function.type = NUTTALL6;
                else if (!strcmp(argument[i], "nuttall7")) fft->window.apodization_function.type = NUTTALL7;
                else if (!strcmp(argument[i], "nuttall8")) fft->window.apodization_function.type = NUTTALL8;
                else if (!strcmp(argument[i], "nuttall9")) fft->window.apodization_function.type = NUTTALL9;
                else if (!strcmp(argument[i], "nuttall10")) fft->window.apodization_function.type = NUTTALL10;
                else if (!strcmp(argument[i], "nuttall11")) fft->window.apodization_function.type = NUTTALL11;
                else if (!strcmp(argument[i], "nuttall12")) fft->window.apodization_function.type = NUTTALL12;
                else if (!strcmp(argument[i], "nuttall14")) fft->window.apodization_function.type = NUTTALL14;
                else if (!strcmp(argument[i], "nuttall15")) fft->window.apodization_function.type = NUTTALL15;
                else if (!strcmp(argument[i], "flat_top")) fft->window.apodization_function.type = FLAT_TOP;
                else if (!strcmp(argument[i], "gaussian")) fft->window.apodization_function.type = GAUSSIAN;
				else if (!strcmp(argument[i], "hanning-poisson")) fft->window.apodization_function.type = HANNING_POISSON;
				else if (!strcmp(argument[i], "helie_a_w1")) fft->window.apodization_function.type = HELIE_A_W1;
				else if (!strcmp(argument[i], "helie_a_w6")) fft->window.apodization_function.type = HELIE_A_W6;
				else if (!strcmp(argument[i], "mod_barlett_hann")) fft->window.apodization_function.type = MOD_BARLETT_HANN;
				else if (!strcmp(argument[i], "bohman")) fft->window.apodization_function.type = BOHMAN;
				else
                {
                    fprintf(stderr, "Unrecognized apodization function: '%s'\n", argument[i]);
                    status = ERROR;
                }
            }
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--klap_ic"))
		{
      		if (++i < n_args) KLAPURI_ITERATION_CONTROL = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--print"))
		{
      		if (++i < n_args)
      		{
            	if (!strcmp(argument[i], "samples")) OPERATION_MODE = PRINT_SAMPLES;
            	else if (!strcmp(argument[i], "f0"))
            	{
            		if (++i < n_args)
            	    {
                        OPERATION_MODE = PRINT_INDEPENDENT_F0_ESTIMATES;
                        if (!strcmp(argument[i], "yin")) INDEP_F0_FUN_PTR = yin_f0_estimate;
                        else if (!strcmp(argument[i], "yan")) INDEP_F0_FUN_PTR = yan_f0_estimate;
                        else if (!strcmp(argument[i], "ac")) INDEP_F0_FUN_PTR = ac_f0_estimate;
                        else if (!strcmp(argument[i], "amdf")) INDEP_F0_FUN_PTR = amdf_f0_estimate;
                        else if (!strcmp(argument[i], "wac")) INDEP_F0_FUN_PTR = wac_f0_estimate;
                        else if (!strcmp(argument[i], "cep")) INDEP_F0_FUN_PTR = cepstrum_f0_estimate;
                        else if (!strcmp(argument[i], "aclos")) INDEP_F0_FUN_PTR = aclos_f0_estimate;
                        else
                        {
                        	fprintf(stderr, "Unrecognized independent F0 technique: '%s'\n", argument[i]);
                        	status = ERROR;
                        	OPERATION_MODE = DEFAULT_OPERATION_MODE;
                        }
                    }
              		else
                	{
                     	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                     	status = ERROR;
                    }
            	}
            	else if (!strcmp(argument[i], "spec"))
             	{
                	if (++i < n_args)
            	    {
                        OPERATION_MODE = PRINT_SPEC;
                        spec_type_set = TRUE;
                        if (!strcmp(argument[i], "ph")) CHOSEN_SPECTRUM = PH_SPECTRUM;
                        else if (!strcmp(argument[i], "uph")) CHOSEN_SPECTRUM = UPH_SPECTRUM;
                        else if (!strcmp(argument[i], "mag")) CHOSEN_SPECTRUM = MAG_SPECTRUM;
                        else if (!strcmp(argument[i], "pow")) CHOSEN_SPECTRUM = POW_SPECTRUM;
                        else if (!strcmp(argument[i], "lp")) CHOSEN_SPECTRUM = LP_SPECTRUM;
                        else if (!strcmp(argument[i], "wd")) CHOSEN_SPECTRUM = WD_SPECTRUM;
                        else if (!strcmp(argument[i], "pn")) CHOSEN_SPECTRUM = PN_SPECTRUM;
                        else
                        {
                        	fprintf(stderr, "Unrecognized spectrum type: '%s'\n", argument[i]);
                        	status = ERROR;
                        	OPERATION_MODE = DEFAULT_OPERATION_MODE;
                            spec_type_set = FALSE;
                        }
                        if (++i < n_args)
                        {
                        	if (!strcmp(argument[i], "doubled_linebreaks"))
                        		DOUBLED_LINEBREAKS = TRUE;
                        	else i--;
                        }
                    }
              		else
                	{
                     	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                     	status = ERROR;
                    }
                }
            	else if (!strcmp(argument[i], "loc_max"))
             	{
                	if (++i < n_args)
            	    {
                        OPERATION_MODE = PRINT_SPEC_LOC_MAX;
                        spec_type_set = TRUE;
                        if (!strcmp(argument[i], "ph")) CHOSEN_SPECTRUM = PH_SPECTRUM;
                        else if (!strcmp(argument[i], "uph")) CHOSEN_SPECTRUM = UPH_SPECTRUM;
                        else if (!strcmp(argument[i], "mag")) CHOSEN_SPECTRUM = MAG_SPECTRUM;
                        else if (!strcmp(argument[i], "pow")) CHOSEN_SPECTRUM = POW_SPECTRUM;
                        else if (!strcmp(argument[i], "lp")) CHOSEN_SPECTRUM = LP_SPECTRUM;
                        else if (!strcmp(argument[i], "wd")) CHOSEN_SPECTRUM = WD_SPECTRUM;
                        else if (!strcmp(argument[i], "pn")) CHOSEN_SPECTRUM = PN_SPECTRUM;
                        else
                        {
                        	fprintf(stderr, "Unrecognized spectrum type: '%s'\n", argument[i]);
                        	status = ERROR;
                        	OPERATION_MODE = DEFAULT_OPERATION_MODE;
                            spec_type_set = FALSE;
                        }
                    }
              		else
                	{
                     	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                     	status = ERROR;
                    }
                }
             	else if (!strcmp(argument[i], "unpred")) OPERATION_MODE = PRINT_UNPRED;
             	else if (!strcmp(argument[i], "threshold")) OPERATION_MODE = PRINT_THRESHOLD;
             	else if (!strcmp(argument[i], "bands_response")) OPERATION_MODE = PRINT_BANDS_RESPONSE;
             	else if (!strcmp(argument[i], "bands_response_summary")) OPERATION_MODE = PRINT_BANDS_RESPONSE_SUMMARY;
             	else if (!strcmp(argument[i], "bands_sum")) OPERATION_MODE = PRINT_BANDS_SUM;
             	else if (!strcmp(argument[i], "onset")) OPERATION_MODE = PRINT_ONSET;
             	else if (!strcmp(argument[i], "f0_estimate")) OPERATION_MODE = PRINT_F0_ESTIMATE;
             	else
              	{
                   	fprintf(stderr, "Unrecognized print mode: '%s'\n", argument[i]);
                   	status = ERROR;
                }
            }
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
        }
        else if (!strcmp(argument[i], "--ruminate"))
        {
        	OPERATION_MODE = RUMINATE;
        }
        else if (!strcmp(argument[i], "--transcribe"))
        {
      		OPERATION_MODE = TRANSCRIBE;
        }
        else if(!strcmp(argument[i], "--piano_roll"))
        {
        	PRINT_NOTE_AS_GNUPLOT_ARROWS = TRUE;
        }
		else
  		{
        	fprintf(stderr, "Unrecognized parameter: '%s'\n", argument[i]);
        	status = ERROR;
        }
	}

	/* validate data */

	if (status == ERROR); /* TODO: verify this if-else block, specially THIS odd line */
	else if (enough_args < 2)
	{
		fprintf(stderr, "Error: insufficient arguments.\n");
		status = ERROR;
	}
	else
	{
        if (!is_power_of_2(fft->window.length)
        	|| fft->window.length < MIN_FFT_WINDOW_SIZE
         	|| fft->window.length > MAX_FFT_WINDOW_SIZE)
        {
            fprintf(stderr, "Error: FFT window size must be a power of 2 between %d "
           			"and %d.\n", MIN_FFT_WINDOW_SIZE, MAX_FFT_WINDOW_SIZE);
            status = ERROR;
        }

        if (!step_set)
        {
        	if (fft->window.apodization_function.type == RECTANGULAR) fft->step = fft->window.length;
        	else fft->step = (unsigned int)(DEFAULT_RELATIVE_FFT_STEP * fft->window.length + 0.5);
        }

        if (fft->zero_padding_ratio > MAX_ZERO_PADDING_RATIO
         	|| !is_power_of_2(fft->zero_padding_ratio))
        {
        	fprintf(stderr, "Error: zero-padding ratio must be a power of two not bigger "
         			"than %d.\n", MAX_ZERO_PADDING_RATIO);
        	status = ERROR;
        }
        if (fft->window.length * fft->zero_padding_ratio > MAX_FFT_WINDOW_SIZE)
        {
        	fprintf(stderr, "Error: total FFT size, i.e. fft_win_len * zero_pad_ratio,\n"
        					"       must be not bigger than %d.\n", MAX_FFT_WINDOW_SIZE);
        	status = ERROR;
        }
        if (fft->step == 0 || fft->step > fft->window.length)
        {
            fprintf(stderr, "Error: window step must be a positive integer not greater "
            	  			"than the FFT window size (%u in this case).\n", fft->window.length);
            /* Note: this is a semantical restriction, not a computational one.
               If the step is bigger than the window size, 'step - size' samples will be
               skipped (i.e. completely ignored) every 'size' samples. */
            status = ERROR;
        }
        if (!verify_extension(input->file_name, wav_ext))
        {
        	fprintf(stderr, "Error: '%s' seems not to be a RIFF WAVE file.\n",
         																	input->file_name);
        	status = ERROR;
        }
        switch (UNPRED_METHOD)
        {
            case LP_UM:
                UNPRED_FUN_PTR = log_power_unpredict;
                break;
            case ILP_UM:
                UNPRED_FUN_PTR = incr_log_power_unpredict;
                break;
            case PH_UM:
                UNPRED_FUN_PTR = phase_unpredict;
                break;
           	case MPH_UM:
	           	UNPRED_FUN_PTR = modified_phase_unpredict;
           		break;
            case COMPLEX_UM:
                UNPRED_FUN_PTR = complex_unpredict;
                break;
            case NON_DECREASING_COMPLEX_UM:
                UNPRED_FUN_PTR = non_decreasing_complex_unpredict;
                break;
            case NEW_UM:
                UNPRED_FUN_PTR = new_unpredict;
                break;
            case RMS_UM:
            	UNPRED_FUN_PTR = rms_unpredict;
            	break;
           	case WRMS_UM:
	           	UNPRED_FUN_PTR = weighted_rms_unpredict;
           		break;
           	case ERB_UM:
	           	UNPRED_FUN_PTR = erb_unpredict;
           		break;
           	case BARK_UM:
	           	UNPRED_FUN_PTR = bark_unpredict;
           		break;
           	case MKL_UM:
	           	UNPRED_FUN_PTR = mkl_unpredict;
           		break;
           	case SPEC_DIF_UM:
	           	UNPRED_FUN_PTR = spectral_difference_unpredict;
           		break;
        }
        switch (ESTIMATION_METHOD)
        {
            case MAX_INDEX:
                ESTIMATION_FUN_PTR = max_index_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_MAX_INDEX_SPECTRUM;
                break;
            case HPS:
                ESTIMATION_FUN_PTR = hps_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_HPS_SPECTRUM;
                break;
            case HSC:
                ESTIMATION_FUN_PTR = hsc_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_HSC_SPECTRUM;
                break;
            case BANDWISE_HSC:
                ESTIMATION_FUN_PTR = bandwise_hsc_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_BANDWISE_HSC_SPECTRUM;
                break;
            case FFT_FFT:
                ESTIMATION_FUN_PTR = fft_fft_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_FFT_FFT_SPECTRUM;
                break;
            case KLAPURI:
                ESTIMATION_FUN_PTR = klapuri_pitch_estimate;
                if (!spec_type_set) CHOSEN_SPECTRUM = DEFAULT_KLAPURI_SPECTRUM;
                break;
        }
    }

    if (status == ERROR)
    {
        fprintf(stderr, "\n"
                "Expected Syntax:\n"
                "  asymut <FFT_SIZE> <ZERO_PADDING_RATIO> <WINDOW_STEP> <file.wav>\n"
                "\n");
    }
    else
    {
    	init_window(&fft->window); /* TODO: move this from here to main() */
    	fft->num_of_bins = fft->window.length * fft->zero_padding_ratio / 2;
    }

    return status;
}


/********************************************************************************************/
/* IS_POWER_OF_2                                                                            */
/********************************************************************************************/
unsigned char is_power_of_2(unsigned long n)
{
    if (!n) return FALSE;

    while (n > 2)
    {
        if (n%2) return FALSE;
        n /= 2;
    }
    return TRUE;
}

/********************************************************************************************/
/* VERIFY_EXTENSION (case insensitive)                                                      */
/********************************************************************************************/
unsigned char verify_extension (const char* name, const char* extension)
{
    unsigned char ext_len = strlen(extension), i;
    unsigned short name_len = strlen(name);

    for (i = 1; i <= ext_len; i++)
    	if (tolower(name[name_len - i]) != tolower(extension[ext_len - i])) return FALSE;

    return TRUE;
}

/********************************************************************************************/
/* OPEN_AUDIOFILE                                                                           */
/********************************************************************************************/
signed char open_audiofile(Audio_Stream* stream)
{
    if ((stream->file_pointer = fopen(stream->file_name, "rb")) == NULL)
    {
        fprintf(stderr, "Error opening '%s'.\n", stream->file_name);
        return ERROR;
    }

    return OK;
}

/********************************************************************************************/
/* CLOSE_AUDIOFILE                                                                          */
/********************************************************************************************/
signed char close_audiofile(Audio_Stream* stream)
{
    if (fclose(stream->file_pointer) == EOF)
        return ERROR;

    return OK;
}

/********************************************************************************************/
/* READ_WAV_HEADER                                                                          */
/********************************************************************************************/
short read_wav_header(Audio_Stream *stream)
{
    long unsigned l_temp = 0, bytes_per_second = 0, bytes_per_sample = 0;
    int unsigned i_temp = 0;
    char s_temp[5], s_ref[5];
    unsigned char t_sample[32];

    s_temp[4] = s_ref[4] = '\0';

    /* block identifier */
    fread(&s_temp, 4, 1, stream->file_pointer);
    strcpy (s_ref, "RIFF");
    if (strcmp(s_temp, s_ref))
    {
        fprintf(stderr, "Error reading input WAV header (RIFF).\n");
        return ERROR;
    }

    /* block size: size of header block (36) + size of samples block
        (number of samples * number of channels * byte depth) */
    l_temp = 0;
    fread(&l_temp, 4, 1, stream->file_pointer);

    /* block identifier */
    fread(&s_temp, 4, 1, stream->file_pointer);
    strcpy (s_ref, "WAVE");
    if (strcmp(s_temp, s_ref))
    {
        fprintf(stderr, "Error reading input WAV header (WAVE).\n");
        return ERROR;
    }

    /* block identifier */
    fread(&s_temp, 4, 1, stream->file_pointer);
    strcpy (s_ref, "fmt ");
    if (strcmp(s_temp, s_ref))
    {
        fprintf(stderr, "Error reading input WAV header (fmt ).\n");
        return ERROR;
    }

    /* format block size */
    l_temp = 0;
    if (fread(&l_temp, 4, 1, stream->file_pointer) < 1)
    {
        fprintf(stderr, "Error reading format block size (pos %ld)\n",
        														ftell(stream->file_pointer));
    }
    if (l_temp != 16)
    {
        fprintf(stderr, "Unsupported WAV format (format block size %ld).\n", l_temp);
        return ERROR;
    }

    /* sample format (PCM = 1) */
    i_temp = 0;
    fread(&i_temp, 2, 1, stream->file_pointer);
    if (i_temp != 1)
    {
        fprintf(stderr, "Unsupported WAV format (sample format = %d).\n", i_temp);
        return ERROR;
    }

    /* number of channels */
    i_temp = 0;
    fread(&i_temp, 2, 1, stream->file_pointer);
    if (i_temp < MIN_NUMBER_OF_CHANNELS || i_temp > MAX_NUMBER_OF_CHANNELS)
    {
        fprintf(stderr, "Unsupported WAV format (number of channels).\n");
        return ERROR;
    }
    stream->num_of_channels = i_temp;

    /* sample rate */
    l_temp = 0;
    fread(&l_temp, 4, 1, stream->file_pointer);
    if (l_temp < MIN_SAMPLE_RATE || l_temp > MAX_SAMPLE_RATE)
    {
        fprintf(stderr, "Unsupported WAV format (sample rate).\n");
        return ERROR;
    }
    stream->sample_rate = l_temp;

    /* bytes per second (sample rate * number of channels * byte depth) */
    fread(&bytes_per_second, 4, 1, stream->file_pointer);

    /* total number of bytes per sample (number of channels * byte depth) */
    fread(&bytes_per_sample, 2, 1, stream->file_pointer);

    /* bit depth */
    i_temp = 0;
    fread(&i_temp, 2, 1, stream->file_pointer);
    stream->bit_depth = i_temp;
    stream->byte_depth = (i_temp/8);

    if (bytes_per_second != stream->sample_rate * stream->num_of_channels
    																	* stream->byte_depth)
    {
        fprintf(stderr, "Inconsistent WAV header (bytes per second).\n");
        return ERROR;
    }
    if (bytes_per_sample != stream->num_of_channels * stream->byte_depth)
    {
        fprintf(stderr, "Inconsistent WAV header (bytes per sample).\n");
        return ERROR;
    }

    /* block identifier */
    fread(&s_temp, 4, 1, stream->file_pointer);
    strcpy (s_ref, "data");
    if (strcmp(s_temp, s_ref))
    {
        fprintf(stderr, "Error reading input WAV header (data).\n");
        return ERROR;
    }

    /* samples block size (number of samples * number of channels * byte depth) */
    l_temp = 0;
    fread(&l_temp, 4, 1, stream->file_pointer);
    stream->total_num_of_samples = l_temp / (stream->num_of_channels * stream->byte_depth);

    /* ensures that the state 'samples block size' correctly reflects reality */
    fseek(stream->file_pointer, (stream->total_num_of_samples - 1)
    							* (stream->byte_depth * stream->num_of_channels), SEEK_CUR);
    if (fread(t_sample, stream->byte_depth * stream->num_of_channels,
    			1, stream->file_pointer) == 0)
    {
        fprintf(stderr, "Invalid samples block size (file is smaller than stated in "
				        "the RIFF WAVE header).\n");
        return( ERROR );
    }
    if (fread(t_sample, stream->byte_depth * stream->num_of_channels,
    			1, stream->file_pointer) != 0)
    {
        fprintf(stderr, "Invalid samples block size (file is bigger than stated in "
				        "the RIFF WAVE header).\n");
        return( ERROR );
    }

    /* place the stream at the beginning of the sample block */
    fseek(stream->file_pointer, 44L, SEEK_SET);

    return OK;
}

/********************************************************************************************/
/* C8_TO_FLOAT                                                                              */
/********************************************************************************************/
float c8_to_float(const unsigned char sample)
{
	return ((float)sample - 128);
}

/********************************************************************************************/
/* C16_TO_FLOAT                                                                             */
/********************************************************************************************/
float c16_to_float(const unsigned char sample[2])
{
    if (sample[1] & 128)
        return (float)(((sample[0] + (sample[1] << 8)) - 1) ^ 0xFFFF) * (-1.0);
    else
        return (float)(sample[0] + (sample[1] << 8));
}

/********************************************************************************************/
/* ROUND_DOUBLE_TO_USHORT                                                                   */
/********************************************************************************************/
unsigned short round_double_to_ushort(double val)
{
    double short_int_val;

    if (val < 0)
    	return 0;

    if (val > USHRT_MAX)
    	return USHRT_MAX;

    return short_int_val + (modf(val, &short_int_val) >= 0.5 ? 1 : 0);
}

/********************************************************************************************/
/* EMALLOC                                                                                  */
/********************************************************************************************/
void* emalloc(const size_t size)
{
    void* p;

    if ( ( p = malloc(size) ) != NULL )
    	return p;

    fprintf(stderr, "Error allocating memory. Aborting...\n");
    exit ( ERROR );
}

/********************************************************************************************/
/* ECALLOC                                                                                  */
/********************************************************************************************/
void* ecalloc(const size_t size)
{
    void* p;

    if ( ( p = calloc(1, size) ) != NULL )
    	return p;

    fprintf(stderr, "Error allocating memory. Aborting...\n");
    exit ( ERROR );
}

/********************************************************************************************/
/* PRINT_FFT_PROP_STATE                                                                     */
/********************************************************************************************/
void print_fft_prop_state(FFT_properties* fft_p, const int n)
{
    fprintf(stderr, "***********************\n");
    fprintf(stderr, "****** PFPS %d ********\n", n);
    fprintf(stderr, "***********************\n");
    fprintf(stderr, "freq_resolution: %g\n",  fft_p->freq_resolution);
    fprintf(stderr, "num_of_bands: %u\n",  fft_p->num_of_bins);
    fprintf(stderr, "step: %u\n",  fft_p->step);
    fprintf(stderr, "time_resolution: %g\n",  fft_p->time_resolution);
    fprintf(stderr, "window.length: %u\n",  fft_p->window.length);
    fprintf(stderr, "zero_padding_ratio: %u\n",  fft_p->zero_padding_ratio);
}

/********************************************************************************************/
/* PRINC_DET                                                                                */
/********************************************************************************************/
double princ_det(double angle)
{
    return fmod(angle, PI);
}

/********************************************************************************************/
/* SELECTION_SORT                                                                           */
/********************************************************************************************/
void selection_sort(double v[], unsigned int n)
{
    int i, aux, max;

    if (n < 2)
        return;

    do
    {
        for (i = 0, max = 0; i < n; i++)
        {
            if (v[i] > v[max])
                max = i;
        }
        if (max != (--n))
        {
            aux = v[n];
            v[n] = v[max];
            v[max] = aux;
        }
    }
    while (n);
}

/********************************************************************************************/
/* CALC_MOV_MEDIAN                                                                          */
/********************************************************************************************/
void calc_mov_median(double signal[], unsigned int sig_len, unsigned short win_len,
							unsigned short med_index, double scaling_factor,
       						double const_part, double median[])
{
    unsigned int i, last_index;
	double window[win_len];

    if (sig_len < win_len)
    {
        fprintf(stderr, "Error: window length must not be greater than signal length.\n");
    	return;
    }
    if (med_index >= win_len)
    {
        fprintf(stderr, "Error: median index must not be greater than window length.\n");
    	return;
   	}

   	last_index = sig_len - win_len;
   	for (i = 0; i <= last_index; i++)
   	{
        memmove(window, signal + i, win_len * sizeof(double));
        selection_sort(window, win_len);
        median[i] = window[med_index] * scaling_factor + const_part;
   	}
}

/********************************************************************************************/
/* PRINT_CURRENT_CONFIG                                                                     */
/********************************************************************************************/
/* TODO: improve the name of this function */
/* TODO: finish this thing */
void print_current_config(Buffer *buf)
{
	printf("# input file: %s\n", buf->stream->file_name);
	printf("# bit depth: %u\n", buf->stream->bit_depth);
	printf("# sample rate: %u\n", buf->stream->sample_rate);
	printf("# number of channels: %u\n", buf->stream->num_of_channels);
	printf("# total number of samples: %lu\n", buf->stream->total_num_of_samples);
	printf("#\n");
	printf("# FFT window length (in samples): %u\n", buf->fft_p->window.length);
	printf("# FFT window step (a.k.a. hop or stride): %u\n", buf->fft_p->step);
	printf("# FFT window zero padding ratio: %u\n", buf->fft_p->zero_padding_ratio);
/*    APODIZATION_FUNCTION;

    OPERATION_MODE;

    UNPRED_METHOD;
    ESTIMATION_METHOD;

	PITCH_RANGE_LOWEST_NOTE;
	PITCH_RANGE_HIGHEST_NOTE;
    MIN_ABSOLUTE_F0;
    MIN_F0_TO_FREQ_RES_RATIO;
    MAX_ABSOLUTE_F0;
    MAX_CRITICAL_BANDS;
    MIN_INTERONSET_DISTANCE;
    ONSET_THRESHOLD_MIN_WINDOW_LENGTH;
    ONSET_THRESHOLD_MAX_WINDOW_LENGTH;
    ONSET_THRESHOLD_PERCENTILE;
    ONSET_THRESHOLD_SCALING_FACTOR;
    ONSET_THRESHOLD_CONST_PART;
    MAX_DELAY_AFTER_ONSET;
    MAX_GAP_INSIDE_NOTE;
    MIN_NOTE_DURATION;

    FREQ_REF_A4;
    POWER_NORMALIZATION_LEVEL;
    KLAPURI_ITERATION_CONTROL;
    */
    printf("#\n");
}

/********************************************************************************************/
/* LOG2                                                                                     */
/********************************************************************************************/
double log2(double x)
{
	return log10(x) / LOGARITHM_OF_2_IN_BASE_10;
}

/********************************************************************************************/
/* HERTZ_TO_BARK                                                                            */
/********************************************************************************************/
/* converts frequency, in Hertz, to Bark Scale band,
 * according to Traunm�ller (1983, 1988, 1990) */
double hertz_to_bark(const double f)
{
	double bark = (26.81*f)/(1960+f)-0.53;

	if (bark < 2) {
		bark = bark + 0.15*(2-bark);
	} else if (bark > 20.1) {
		bark = bark + 0.22*(bark-20.1);
	}

	return bark;
}

/********************************************************************************************/
/* HERTZ_TO_ERB                                                                             */
/********************************************************************************************/
/* converts frequency, in Hertz, to Equivalent Rectangular Bandwidth band */
double hertz_to_erb(const double f)
{
	return 21.4 * log10(4.37e-3*f + 1);
}

/********************************************************************************************/
/* ERB_TO_HERTZ                                                                             */
/********************************************************************************************/
/* converts Equivalent Rectangular Bandwidth band to frequency in Hertz */
double erb_to_hertz(const double b)
{
	return (pow(10, b/21.4) - 1) / 4.37e-3;
}
