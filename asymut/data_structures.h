#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

#include <stdio.h> /* needed because of FILENAME_MAX */
#include "complex.h"

typedef enum {FALSE, TRUE} Boolean;

typedef enum {PRINT_WINDOW_TIME_RESPONSE, PRINT_WINDOW_FREQ_RESPONSE, PRINT_BANDS_RESPONSE,
	PRINT_BANDS_RESPONSE_SUMMARY, PRINT_BANDS_SUM, PRINT_INDEPENDENT_F0_ESTIMATES,
	PRINT_SAMPLES, PRINT_SPEC, RUMINATE, PRINT_SPEC_LOC_MAX, PRINT_UNPRED, PRINT_THRESHOLD,
	PRINT_ONSET, PRINT_F0_ESTIMATE, TRANSCRIBE } Asymut_Operation_Mode;
   	/* Warning: DO NOT change this order */

typedef enum {MAX_INDEX, HPS, HSC, BANDWISE_HSC, FFT_FFT, KLAPURI} F0_Estimation_Method;

typedef enum {LP_UM, ILP_UM, PH_UM, MPH_UM, COMPLEX_UM, NON_DECREASING_COMPLEX_UM, NEW_UM, RMS_UM,
				WRMS_UM, ERB_UM, BARK_UM, MKL_UM, SPEC_DIF_UM} Unpredictability_Method;

typedef enum {PH_SPECTRUM, UPH_SPECTRUM, MAG_SPECTRUM, POW_SPECTRUM, LP_SPECTRUM,
	WD_SPECTRUM, PN_SPECTRUM} Spectrum_Type;

typedef enum {TRIANGULAR, RECTANGULAR, BLACKMAN_HARRIS, FLAT_TOP, HANN, BLACKMAN,
	NUTTALL3, NUTTALL4, NUTTALL5, NUTTALL6, NUTTALL7, NUTTALL8, NUTTALL9, NUTTALL10,
	NUTTALL11, NUTTALL12, HAMMING, NUTTALL14, NUTTALL15, GAUSSIAN, HANNING_POISSON,
	HELIE_A_W1,	HELIE_A_W6, MOD_BARLETT_HANN, BOHMAN} Apodization_Function_Type;


typedef struct t_audio_stream
{
    unsigned int sample_rate;
    unsigned char bit_depth;
    unsigned char byte_depth;
    unsigned char num_of_channels;
    unsigned long total_num_of_samples; /* per channel */
    char file_name[FILENAME_MAX];
    FILE* file_pointer;
}
Audio_Stream;


typedef struct t_Apodization_Function
{
	Apodization_Function_Type type;

	double *parameters;
	unsigned char num_of_parameters;
}
Apodization_Function;


typedef struct t_window
{
	Apodization_Function apodization_function;

    unsigned int length;
    double* content; /* TODO: consider changing the name to 'apodization function' or other */
    double* derivative;

    double total_area;
}
Window;


typedef struct t_FFT_properties
{
    Window window;
    unsigned short zero_padding_ratio;

    unsigned short num_of_bins;
    double freq_resolution; /* 'sample rate/(window.size * zero_padding_ratio)' */
    double time_resolution; /* 'step/sample rate' */

    /* TODO: refactor 'step' to 'hop' (including internal function variables) */
    unsigned short step; /* 'window.size * (1 - overlap)', a.k.a 'hop' or 'stride' */
}
FFT_properties;


typedef struct t_spectrum
{
    double* content;

/*	WARNING: BE CAREFUL NOT TO FORGET CALCULATING AND RECYCLING THE FOLLOWING EXTRA DATA */

/*	double* accumulated; */
/*	double* median; */
    double total_sum;

    unsigned short max_index;
/*	unsigned short min_index; */
}
Spectrum;


typedef struct t_spectra
{
	Spectrum phase;
	Spectrum unwrapped_phase;
    Spectrum magnitude;
    Spectrum power;
    Spectrum log_power;
    Spectrum warped_denoised;
    Spectrum power_noise;
}
Spectra;



/* TODO: replace this struct by something more refined (estudo.c) */
typedef struct t_note
{
	char name[3];
    unsigned char oct;
    unsigned char midi_number;

    double fund_freq;

    unsigned char state; /* ON, OFF */

    double onset_time; /* in seconds */
    double last_evidence_time; /* in seconds */

/*	double last_salience_level; */ /* even if not sounding, note has a small level of presence
 * 									  (something like its probability of being present) */
}
Note;


typedef struct t_pitch_range
{
    unsigned char num_of_notes;
    unsigned char current_polyphony;

    Note* note;
}
Pitch_Range;

typedef struct t_bands_properties
{
	unsigned char num_of_bands;
	double *freq_bound;
}
Bands_properties;

typedef struct t_buffer
{
	FFT_properties* fft_p;
	Audio_Stream* stream;
	unsigned short threshold_win_len;

    unsigned int sample_array_length; 				/* total samples per channel */
    float** sample_data;        						/* multichannel */

    unsigned short spec_to_sample_offset;
    unsigned int spec_array_length;
    Spectra** spec_data;        						/* multichannel */

    unsigned short unpredictability_to_sample_offset;
    unsigned short unpredictability_to_spec_offset;
    unsigned int unpredictability_array_length;
    double* unpredictability_data;   					/* SINGLE channel */

    unsigned int threshold_to_sample_offset;
    unsigned short threshold_to_spec_offset;
    unsigned short threshold_to_unpredictability_offset;
    unsigned int threshold_array_length;
    double* threshold_data;								/* SINGLE channel */

    unsigned int actual_num_of_samples;
    unsigned int actual_num_of_specs;
    unsigned int actual_num_of_unpredicts;
    unsigned int actual_num_of_thresholds;

    unsigned int trailing_zeroes_left;

    signed long first_sample_file_pos;
}
Buffer;


typedef struct t_critical_band
{
    unsigned short first_bin; /* first non null bin */
    unsigned short num_of_bins; /* total number of non null bins */
    double center; /* in Hertz */
    double width; /* in musical octaves */
    double* response;
    double response_total_sum;
}
Critical_Band;


typedef struct t_auditory_model
{
    FFT_properties* fft_p;

    unsigned char num_of_bands;
    Critical_Band* band;
}
Auditory_Model;

typedef void (*Postponable_Function)(void);

#endif
