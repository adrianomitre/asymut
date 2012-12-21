/* ASYMUT (Automated SYstem of Music Transcription) - author: Adriano Brito Mitre */
/* TODO: consider changing the project name */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "data_structures.h"
#include "numeric_constants.h"
#include "config_variables.h"
#include "aux_fun.h"
#include "buffer.h"
#include "spec_fun.h"
#include "notes.h"

#include "postponing.h"

#include "stdft.h"
#include "trigonometric_constants.h"

int main(int argc, char *argv[])
{
	/*
	const unsigned int n = 8, w = 0, z = n;
	Complex *s, *spec;
	double *rs, *pow, *ph;

	rs = (double*)emalloc(n*sizeof(double));
	pow = (double*)emalloc((z+n/2)*sizeof(double));
	ph = (double*)emalloc((z+n/2)*sizeof(double));

	s = (Complex*)emalloc(sizeof(Complex)*n);
	spec = (Complex*)emalloc(sizeof(Complex)*n);
	unsigned int i;

	for (i = 0; i < n; i++) {
		rs[i] = s[i].re = cos(w*TWO_PI/(double)n*i);
		s[i].im = 0;
	}
	reference_dft(s, spec, n);
	for (i = 0; i < n; i++) {
		double mag, ph;
		c_complex_to_polar(spec[i], &mag, &ph);
		printf("%d %g %g\n", i, mag*mag, ph);
	}
	reference_fft_real_signal_power_phase(rs, n, z, pow, ph);
	for (i = 0; i < (z+n)/2; i++) {
		printf("%d %g %g\n", i, pow[i], ph[i]);
	}
	return 0;
	*/

	/* TODO: consider initializing these two with NULL to be allocated inside
	 * parse_arguments() or, immediately before it, with something like init_audio_stream() */
    Audio_Stream* stream = (Audio_Stream*) emalloc(sizeof(Audio_Stream));
    FFT_properties* fft_prop = (FFT_properties*) emalloc(sizeof(FFT_properties));

    Buffer* buf = NULL;
    Auditory_Model* aud_mod = NULL;
    Pitch_Range* pr = NULL;

    /****************************************************************************************/

    prepare_postponing(MAX_POSTPONED_FUNCTIONS);

    if (parse_arguments(argc, argv, stream, fft_prop) != OK)
    	return ERROR;

    if (open_audiofile(stream) != OK)
    	return ERROR;
    if (read_wav_header(stream) != OK)
    	return ERROR;

    if (init_buffer(&buf, stream, fft_prop) != OK)
        return ERROR;

	if (OPERATION_MODE >= PRINT_SPEC) {
		init_normalization(buf);
	    calc_spectrum_offset(buf, CHOSEN_SPECTRUM);
	}

    if (OPERATION_MODE >= PRINT_BANDS_RESPONSE)
    {
        if (init_pitch_range(&pr, PITCH_RANGE_LOWEST_NOTE, PITCH_RANGE_HIGHEST_NOTE) != OK)
        	return ERROR;

        init_auditory_model(&aud_mod, buf);
    }

	if (OPERATION_MODE >= PRINT_SPEC && OPERATION_MODE != RUMINATE) {
		print_current_config(buf);
	}

    /****************************************************************************************/

	switch (OPERATION_MODE) {
		case PRINT_INDEPENDENT_F0_ESTIMATES:
		    MIN_ALLOWED_LAG = ceil(buf->stream->sample_rate/MAX_ABSOLUTE_F0);
			MAX_ALLOWED_LAG = floor(buf->stream->sample_rate/MIN_ABSOLUTE_F0);
			while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
	    		print_buf_indep_f0_estimates(buf, INDEP_F0_FUN_PTR);
			break;
		case PRINT_SAMPLES:
			while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
	    		print_buf_samples(buf);
			break;
		case PRINT_WINDOW_TIME_RESPONSE:
			break;
		case PRINT_WINDOW_FREQ_RESPONSE:
			break;
		case PRINT_BANDS_RESPONSE:
			 print_auditory_model(aud_mod);
			 break;
		case PRINT_BANDS_RESPONSE_SUMMARY:
	        print_auditory_model_summary(aud_mod);
	        break;
	    case PRINT_BANDS_SUM:
			print_auditory_model_per_bin_sum(aud_mod);
			break;
		case PRINT_SPEC:
	    	while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
    			print_buf_spec(buf, print_spec);
			break;
		case PRINT_SPEC_LOC_MAX:
	    	while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
        print_buf_spec(buf, print_spec_interp_grandke);
				/*
				print_buf_spec(buf, print_spec_interp_parabolic);
        */

			break;
		case RUMINATE:
			break;
		case PRINT_UNPRED:
	        while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
    	        print_buf_unpred(buf);
			break;
		case PRINT_THRESHOLD:
	        while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
            	print_buf_threshold(buf);
			break;
		case PRINT_ONSET:
	    	while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
    	        print_buf_onset(buf);
			break;
		case PRINT_F0_ESTIMATE:
	    	while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
    	        print_buf_pitch_estimate(buf, aud_mod, pr, ESTIMATION_FUN_PTR);
			break;
		case TRANSCRIBE:
	    	if (PRINT_NOTE_AS_GNUPLOT_ARROWS) printf("unset arrow\n");
        	while ( fill_buffer(buf, stream, aud_mod, UNPRED_FUN_PTR) == OK )
     	       transcribe_buf(buf, aud_mod, pr, ESTIMATION_FUN_PTR);
			break;
		default:
	        fprintf(stderr, "Error: unknown operation mode\n");
    	    free_buffer(buf);
        	return ERROR;
	}

    /****************************************************************************************/

    free_buffer(buf);

    if (close_audiofile(stream) != OK)
    	return ERROR;

    call_postponed();
    finish_postponing();

    return OK;
}
