#include "data_structures.h"

signed char init_buffer(Buffer **buffer, Audio_Stream* stream, FFT_properties *fft_p);

signed char fill_buffer(Buffer* buf, Audio_Stream* input_stream, Auditory_Model* aud_mod,
								double unpredict_method(Buffer*, unsigned char, unsigned int));

signed char is_buffer_full_of_samples(Buffer *buf);
signed char is_buffer_full_of_thresholds(Buffer *buf);

signed char read_file_to_buffer(Audio_Stream* input_stream, Buffer* buf);
void recycle_buffer(Buffer* buf);
void pad_with_trailing_zeros(Buffer* buf);

signed char update_buffer(Buffer* buf, Auditory_Model* aud_mod,
						double unpredict_method(Buffer*, unsigned char, unsigned int));

void calc_threshold(Buffer* buf, unsigned short med_index, double scale_factor,
																		double const_part);

void transcribe_buf(Buffer* buf, Auditory_Model* aud_mod, Pitch_Range* pr,
				void (*pitch_estimate)(Spectrum*, Auditory_Model*, Pitch_Range*, double));

void free_buffer(Buffer* buf);
signed char is_buffer_consistent(Buffer* buf, FFT_properties* fft_p);
void print_buffer_state(Buffer* buf, const int n);
void print_buf_indep_f0_estimates(Buffer* buf,
					void (*f0_method)(const float*, const unsigned int, double*, double*));
void print_buf_samples(Buffer* buf);

void print_buf_spec(Buffer* buf,
					void (*print_spec)(Buffer*, const unsigned int, const unsigned int));

void print_buf_unpred(Buffer* buf);
void print_buf_threshold(Buffer* buf);
void print_buf_onset(Buffer* buf);

void print_buf_pitch_estimate(Buffer* buf, Auditory_Model* aud_mod, Pitch_Range* pr,
					void (*pitch_estimate)(Spectrum*, Auditory_Model*, Pitch_Range*, double));
