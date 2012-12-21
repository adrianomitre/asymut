#ifndef _CONFIG_CONSTANTS_H_
#define _CONFIG_CONSTANTS_H_

#include "data_structures.h"

double static const TEMPORAL_RESOLUTION = 0.00005; /* in seconds */

Note_Printing_Style static const DEFAULT_NOTE_PRINTING_STYLE = MIDIFIABLE_CSV;

unsigned char static const DEFAULT_GM_INSTRUMENT_NUMBER = 17;
unsigned char static const DEFAULT_GM_CHANNEL_VOLUME = 127;

unsigned char static const DEFAULT_PITCH_RANGE_LOWEST_NOTE = 20; /* in MIDI note number */
unsigned char static const DEFAULT_PITCH_RANGE_HIGHEST_NOTE = 109; /* in MIDI note number */

double static const DEFAULT_MAX_DELAY_AFTER_ONSET = 0.046; /* in seconds */
double static const DEFAULT_MAX_GAP_INSIDE_NOTE = 0.036; /* in seconds */
double static const DEFAULT_MIN_NOTE_DURATION = 0.05; /* in seconds */

double static const DEFAULT_MAX_ATTACK_LENGTH = 0.2; /* in seconds */
double static const DEFAULT_MIN_ATTACK_ABSOLUTE_INCREASE = 200;
double static const DEFAULT_MIN_ATTACK_RELATIVE_INCREASE = 0.2;

double static const DEFAULT_FREQ_REF_A4 = 440; /* in Hertz */

double static const DEFAULT_TEMPO_IN_BEATS_PER_MINUTE = 120;
double static const DEFAULT_INITIAL_OFFSET_IN_SECONDS = 0;

double static const DEFAULT_MAX_INTENSITY_LEVEL = 1200; /* setting to 0 disables */
double static const DEFAULT_MIN_RELATIVE_INTENSITY_LEVEL = 0.0; /* setting to 0 disables */
double static const DEFAULT_MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL = 0.0; /* setting to 0 disables */

unsigned int static const MAX_LINE_LENGTH = 1024; /* in characters */

char static const COMMENT_CHARACTER = '#';

Boolean static const DEFAULT_NOTES_WITHOUT_ONSET_PERMISSION = TRUE;

#endif /* _CONFIG_CONSTANTS_H_ */
