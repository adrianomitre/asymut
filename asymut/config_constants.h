#ifndef _CONFIG_CONSTANTS_H_
#define _CONFIG_CONSTANTS_H_

#include <float.h>
#include "data_structures.h"

enum {MIN_NUMBER_OF_CHANNELS = 1, MAX_NUMBER_OF_CHANNELS = 2};
enum {MIN_BIT_DEPTH = 8, MAX_BIT_DEPTH = 16};
enum {MIN_SAMPLE_RATE = 8000, MAX_SAMPLE_RATE = 48000};
enum {MAX_TOTAL_FFT_SIZE = 65536}; /* TOTAL_FFT = FFT_WIN_SIZE * ZERO_PADDING_RATIO */
enum {MIN_FFT_WINDOW_SIZE = 32, MAX_FFT_WINDOW_SIZE = 65536};
enum {MAX_ZERO_PADDING_RATIO = MAX_TOTAL_FFT_SIZE/MIN_FFT_WINDOW_SIZE};

unsigned char static const DEFAULT_PITCH_RANGE_LOWEST_NOTE = 21; /* in MIDI note number */
unsigned char static const DEFAULT_PITCH_RANGE_HIGHEST_NOTE = 108; /* in MIDI note number */

double static const DEFAULT_MIN_ABSOLUTE_F0 = 40; /* in Hertz */
double static const DEFAULT_MIN_F0_TO_FREQ_RES_RATIO = 0; /* zero-padding is not considered
						 setting to 0 disables, common values: 1.9 (local maxima harmonic detection),
              										 						2.2 ~ 2.6 (harmonic interpolation without zero-padding),
              										 						3.3, 4.5 */
double static const DEFAULT_MAX_ABSOLUTE_F0 = 4000; /* in Hertz */ /* previously 4308.668 */
unsigned char static const DEFAULT_MAX_CRITICAL_BANDS = 18;

double static const DEFAULT_MIN_INTERONSET_DISTANCE = 0.05; /* in seconds */

double static const DEFAULT_ONSET_THRESHOLD_MIN_WINDOW_LENGTH = 0.20; /* in seconds */
double static const DEFAULT_ONSET_THRESHOLD_MAX_WINDOW_LENGTH = 0.30; /* in seconds */
double static const DEFAULT_ONSET_THRESHOLD_PERCENTILE = 0.5; /* in per-unit, rounded */
double static const DEFAULT_ONSET_THRESHOLD_SCALING_FACTOR = 1.15;
double static const DEFAULT_ONSET_THRESHOLD_CONST_PART = 4.5e3;

double static const DEFAULT_MAX_DELAY_AFTER_ONSET = 0.08; /* in seconds */
double static const DEFAULT_MAX_GAP_INSIDE_NOTE = 0.02; /* in seconds */
double static const DEFAULT_MIN_NOTE_DURATION = 0.1; /* in seconds */

double static const DEFAULT_FREQ_REF_A4 = 440; /* in Hertz */

unsigned char static const DEFAULT_ZERO_PADDING_RATIO = 1;
										
double static const DEFAULT_RELATIVE_FFT_STEP = 0.5; /* in window fraction, i.e.,pertains to (0,1] */
															/* exceptionally, the rectangular window has no overlapping by default */
															/* previous value: 0.390625 */
double static const DEFAULT_NORMALIZATION_LEVEL = 96; /* in decibels */
Boolean static const DEFAULT_NORMALIZATION_ENABLED = TRUE; /* in decibels */
Apodization_Function_Type static const DEFAULT_APODIZATION_FUNCTION_TYPE = HANN;

double static const DEFAULT_KLAPURI_ITERATION_CONTROL = 0.8;

Unpredictability_Method static const DEFAULT_UNPRED_METHOD = ERB_UM;
F0_Estimation_Method static const DEFAULT_ESTIMATION_METHOD = KLAPURI;
Asymut_Operation_Mode static const DEFAULT_OPERATION_MODE = PRINT_SPEC_LOC_MAX;

Spectrum_Type static const DEFAULT_MAX_INDEX_SPECTRUM = MAG_SPECTRUM;
Spectrum_Type static const DEFAULT_HPS_SPECTRUM = MAG_SPECTRUM;
Spectrum_Type static const DEFAULT_HSC_SPECTRUM = LP_SPECTRUM;
Spectrum_Type static const DEFAULT_BANDWISE_HSC_SPECTRUM = MAG_SPECTRUM;
Spectrum_Type static const DEFAULT_FFT_FFT_SPECTRUM = MAG_SPECTRUM;
Spectrum_Type static const DEFAULT_KLAPURI_SPECTRUM = WD_SPECTRUM;

double static const DEFAULT_GAUSSIAN_WINDOW_ALPHA = 3.0;
double static const DEFAULT_HANNING_POISSON_WINDOW_ALPHA = 2.0;

double static const MIN_AUDIBLE_FREQUENCY = 20;
double static const MAX_AUDIBLE_FREQUENCY = 20e3;

unsigned int static const MAX_POSTPONED_FUNCTIONS = 256;

#endif /* _CONFIG_CONSTANTS_H_ */
