#ifndef _CONFIG_VARIABLES_H_
#define _CONFIG_VARIABLES_H_

#include "config_constants.h"
#include "data_structures.h"

unsigned char PITCH_RANGE_LOWEST_NOTE; /* in MIDI note number */
unsigned char PITCH_RANGE_HIGHEST_NOTE; /* in MIDI note number */

double MIN_ABSOLUTE_F0; /* in Hertz */
double MIN_F0_TO_FREQ_RES_RATIO; /* zero-padding is not considered
													 setting to 0 disables
              										 common values: 2.4, 3.3 */
double MAX_ABSOLUTE_F0; /* in Hertz */
unsigned char MAX_CRITICAL_BANDS;

double MIN_INTERONSET_DISTANCE; /* in seconds */

double ONSET_THRESHOLD_MIN_WINDOW_LENGTH; /* in seconds */
double ONSET_THRESHOLD_MAX_WINDOW_LENGTH; /* in seconds */
double ONSET_THRESHOLD_PERCENTILE; /* in per-unit (i.e., percent * 100)
														rounde to the  */
double ONSET_THRESHOLD_SCALING_FACTOR;
double ONSET_THRESHOLD_CONST_PART;

double MAX_DELAY_AFTER_ONSET; /* in seconds */
double MAX_GAP_INSIDE_NOTE; /* in seconds */

double MIN_NOTE_DURATION; /* in seconds */

double FREQ_REF_A4; /* in Hertz */

double KLAPURI_ITERATION_CONTROL;

Asymut_Operation_Mode OPERATION_MODE;

F0_Estimation_Method ESTIMATION_METHOD;

Unpredictability_Method UNPRED_METHOD;

double (*UNPRED_FUN_PTR)(Buffer*, unsigned char, unsigned int);

void (*ESTIMATION_FUN_PTR)(Spectrum*, Auditory_Model*, Pitch_Range*, double);

void (*INDEP_F0_FUN_PTR)(const float*, const unsigned int, double*, double*);

double NORMALIZATION_LEVEL; /* in decibels */

Boolean NORMALIZATION_ENABLED;

Spectrum_Type CHOSEN_SPECTRUM;

unsigned int CHOSEN_SPECTRUM_OFFSET;

unsigned char DOUBLED_LINEBREAKS;

unsigned char PRINT_NOTE_AS_GNUPLOT_ARROWS;

unsigned int MIN_ALLOWED_LAG;
unsigned int MAX_ALLOWED_LAG;

#endif

/*
#define EBUG
*/

#ifdef EBUG
#define DEBUG(A) A
#else
#define DEBUG(A)
#endif
