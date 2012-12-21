#ifndef _STDFT_H_
#define _STDFT_H_

#include "data_structures.h"
#include "complex.h"

void reference_dft(const Complex s[], Complex ft[], const unsigned long n);
unsigned long unsigned_index(const signed long n, const unsigned long window_size);

void fft(const Complex in[], Complex out[], const unsigned int n);

void fft_real_signal_to_mag(double signal[], FFT_properties* fft_properties,
																		double magnitude[]);

void reference_fft_real_signal_power_phase(const float s[], const unsigned int m,
				const Apodization_Function_Type apod, const unsigned int zeroes,
				double pow[], double ph[]);

void old_prepare_fft_real_signal_to_mag_phase(FFT_properties* fft_properties);
void old_finalize_fft_real_signal_to_mag_phase(void);

double old_fft_real_signal_to_mag_phase(double signal[], FFT_properties* fft_properties,
														double magnitude[], double phase[]);


void prepare_fft_real_signal_to_power_phase(const unsigned int n);
void finalize_fft_real_signal_to_power_phase(void);

void prepare_fft_real_signal(const unsigned int n);
void finalize_fft_real_signal(void);

void fft_real_signal_to_power_phase(double power[], double phase[],
									const double signal[], const unsigned int n);

void fft_real_signal(Complex spectrum[], const double signal[], const unsigned int n);

void init_rev_table(const unsigned int n);
void free_rev_table(void);
void init_sin_table(const unsigned int n);
void free_sin_table(void);
void init_cos_table(const unsigned int n);
void free_cos_table(void);

void init_window(Window *w);
void print_window_time_response(Window* w);
void print_window_freq_response(FFT_properties *fft_p, unsigned int max_bin);
void free_window(Window* w);

void triangular(Window* w);

void rectangular(Window* w);
void blackman_harris(Window* w);
void flat_top(Window* w);
void hann(Window* w);
void blackman(Window* w);
void nuttall3(Window* w); /* a.k.a. "exact" Blackman */
void nuttall4(Window* w);
void nuttall5(Window* w);
void nuttall6(Window* w);
void nuttall7(Window* w);
void nuttall8(Window* w);
void nuttall9(Window* w);
void nuttall10(Window* w);
void nuttall11(Window* w);
void nuttall12(Window* w);
void hamming(Window* w);
void nuttall14(Window* w);
void nuttall15(Window* w);

void mod_barlett_hann(Window* w);
void bohman(Window* w);

void raised_cosine(Window* w, const double coef[], const unsigned char num_of_coef);

void gaussian(Window* w);
void gaussian_param(Window* w, double alpha);

void hanning_poisson(Window* w);
void hanning_poisson_param(Window* w, double alpha);

void helie_a_w1(Window* w);
void helie_a_w6(Window* w);
void helie_a_param(Window* w, double a, double b);

double raised_cosine_mag_resp(const double f, Window *w);

double hann_win_log_power_freq_resp(double f);
double hann_win_mag_freq_resp(double f);
double rectangular_win_mag_freq_resp(double f);

void hann_derivative(Window* w);

double window_mag_resp(const double freq, Window* window);
double window_log_power_resp(const double freq, Window* window);

/* no longer being used, for the sake of computational speed */
void reference_fft(const Complex in[], Complex out[], const unsigned int n);
unsigned int rev(unsigned int k, unsigned int n);

#endif /* _STDFT_H_ */
