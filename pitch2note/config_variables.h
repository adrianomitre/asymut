#ifndef _CONFIG_VARIABLES_H_
#define _CONFIG_VARIABLES_H_

#include "data_structures.h"

unsigned char PITCH_RANGE_LOWEST_NOTE; /* in MIDI note number */
unsigned char PITCH_RANGE_HIGHEST_NOTE; /* in MIDI note number */

double MAX_DELAY_AFTER_ONSET; /* in seconds */
double MAX_GAP_INSIDE_NOTE; /* in seconds */

double MIN_NOTE_DURATION; /* in seconds */

double MAX_ATTACK_LENGTH; /* Maximum time it takes from the first note */
	/* detection until it reaches its full power. After this time, overcoming */
	/* the max intensity will produce a new note. */
double MIN_ATTACK_ABSOLUTE_INCREASE;
double MIN_ATTACK_RELATIVE_INCREASE;
	

double FREQ_REF_A4; /* in Hertz */

double MAX_INTENSITY_LEVEL; /* setting to 0 disables */
double MIN_RELATIVE_INTENSITY_LEVEL; /* setting to 0 disables */
double MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL; /* setting to 0 disables */

Note_Printing_Style PRINT_NOTE_AS;

double TEMPO; /* in beats-per-minute */
double INITIAL_PAUSE; /* in seconds */
unsigned char INSTRUMENT; /* General-MIDI instrument number */
unsigned char CHANNEL_VOLUME; /* General-MIDI channel volume */

Boolean ALLOW_NOTES_WITHOUT_ONSET;

#endif /* _CONFIG_VARIABLES_H_ */

/*
#define EBUG
*/

#ifdef EBUG
#define DEBUG(A) A
#else
#define DEBUG(A)
#endif
