#ifndef _CONFIG_CONSTANTS_H_
#define _CONFIG_CONSTANTS_H_

char static const COMMENT_CHARACTER = '#';

const static double MAX_PEAKS = 32768;
/*
const unsigned int MAX_STRING = 256;
*/

const static double DEFAULT_MIN_F0 = 26.717; /* in Hertz */
const static double DEFAULT_MAX_F0 = 5000; /* in Hertz */ /* previously 4308.668 */

unsigned int static const MAX_LINE_LENGTH = 1024; /* in characters */

const static double MIN_HARM_TO_FUND_RATIO = 0.85;
const static double MAX_PERCENTILE_FREQ_CANDIDATE = 1.00; /* to prevent higher pitch errors */
const static double MIN_CANDIDATE_TO_TOTAL_RATIO = 0.00;
const static double MIN_SUM_RATIO = 0.5;
const static double MIN_ODD_AVERAGE2_RATIO = 0.66;
const static double MIN_AVERAGE2_RATIO = 0.90;

const static double DEFAULT_MIN_RELATIVE_ENERGY = 1; /* default 0.9, older 47.0/60.0; */
const static double DEFAULT_MAX_RELATIVE_ENERGY_ABOVE = 1;
const static double DEFAULT_MAX_RELATIVE_ENERGY_BELOW = 1;

const static unsigned int MAX_F0_CANDIDATES = 1024;

#endif /* _CONFIG_CONSTANTS_H_ */