#ifndef _ONSET_DETECTION_H_
#define _ONSET_DETECTION_H_

#include "data_structures.h"

double log_power_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double incr_log_power_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double phase_unpredict(Buffer* buf, unsigned char ch, unsigned int index);
double modified_phase_unpredict(Buffer* buf, unsigned char ch, unsigned int index);
double complex_unpredict(Buffer* buf, unsigned char ch, unsigned int index);
double non_decreasing_complex_unpredict(Buffer* buf, unsigned char ch, unsigned int index);
double new_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double rms_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double erb_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double bark_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double mkl_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double weighted_rms_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);
double spectral_difference_unpredict(Buffer* buf, unsigned char chan, unsigned int unpred_index);


#endif /* _ONSET_DETECTION_H_ */
