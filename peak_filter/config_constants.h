#ifndef _CONFIG_CONSTANTS_H_
#define _CONFIG_CONSTANTS_H_

#include <float.h>

const static double MAX_NUM_OF_PEAKS = 16384;
const static unsigned int MAX_STRING_LENGTH = 256;

const static double DEFAULT_MIN_FREQUENCY = 20;
const static double DEFAULT_MAX_FREQUENCY = 20e3;
const static double DEFAULT_MIN_ABSOLUTE_MAGNITUDE = 0; /* formerly DBL_MIN */
const static double DEFAULT_MAX_MAGNITUDE_GAP = DBL_MAX; /* formerly 40 */
const static double DEFAULT_MAX_ADJACENT_MAGNITUDE_GAP = DBL_MAX; /* formerly 9 */

char static COMMENT_CHARACTER = '#';
unsigned int static const MAX_LINE_LENGTH = 1024; /* in characters */

/* out of use */
const static double DEFAULT_LOW_MASKING_LIMIT = 1; /* 2.08; */
const static double DEFAULT_HIGH_MASKING_LIMIT = 1; /* 4.82; */


#endif /* _CONFIG_CONSTANTS_H_ */
