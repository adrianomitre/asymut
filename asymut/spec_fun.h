#ifndef _SPEC_FUN_H_
#define _SPEC_FUN_H_

#include "data_structures.h"

void yin_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);
void yan_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);
void ac_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);						
void amdf_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);
void wac_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);
void cepstrum_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);
void aclos_f0_estimate(const float s[], const unsigned int frame_len,
						double *period_estimate, double *estimate_quality);


void copy_spectra(Buffer* buf, const unsigned char chan, const unsigned int src,
																	const unsigned int dest);

void copy_spectrum(Spectrum *src, Spectrum *dest, const unsigned int num_of_bins);

void fft_ph_uph_mag_pow_lp(Buffer* buf, const unsigned char channel,
						const unsigned int sample_index, const unsigned int spec_index);

void fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod);
	
void new_fft_ph_uph_mag_pow_lp_wd_pn(Buffer* buf, const unsigned char channel,
	const unsigned int sample_index, const unsigned int spec_index, Auditory_Model* aud_mod);
	
double mag_warp(double original[], double warped[], const unsigned int min,
											const unsigned int max, const unsigned int n);

void noise_estimate_and_supress(double spec[], double power_noise[], double warp_factor,
																	Auditory_Model* aud_mod);

void mag_warp_denoise(double power_spec[], double warped_denoised[], double power_noise[],
																	Auditory_Model* aud_mod);

void print_spec(Buffer* buf, const unsigned int spec_index, unsigned int offset);

void print_spec_interp(Buffer* buf, const unsigned int spec_index, unsigned int offset);

void print_unpred(Buffer* buf, unsigned int index);
void print_threshold(Buffer* buf, unsigned int index);

void init_auditory_model(Auditory_Model **auditory_model, Buffer* buf);
void print_auditory_model(Auditory_Model* aud_mod);
void print_auditory_model_summary(Auditory_Model* aud_mod);
void print_auditory_model_per_bin_sum(Auditory_Model* aud_mod);
void print_auditory_model_per_band_sum(Auditory_Model* aud_mod);
void print_auditory_model_per_band_average(Auditory_Model* aud_mod);
void print_auditory_model_limits(Auditory_Model* aud_mod);

void klapuri_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void max_index_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void hps_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void hsc_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void hsc_pitch_estimate_2(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void hsc_pitch_estimate_3(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void bandwise_hsc_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void fft_fft_pitch_estimate(Spectrum* spec, Auditory_Model* aud_mod, Pitch_Range* pr,
																			double curr_time);

void calc_spectrum_offset(Buffer *buf, Spectrum_Type);

void calc_spectra_extra_data(Spectra* spec, unsigned int num_of_bins);

double grandke_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);

double barycentric_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);

double parabolic_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);
														
double amplitude_parabolic_interp(double *log_power_spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);
														
double centroid_interp(double *v, const unsigned int n, const unsigned int center);

double quinn_2nd_interp(Spectra *spectra, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);

double quinn_2nd_interp_polynomial(const double x);

void print_spec_interp_quinn(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);

double derivative_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp);
						
double reassignment_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp);

double agrez_interp(double *spec_content, const unsigned int num_of_bins,
														const unsigned int bin_to_interp);

void print_spec_interp_parabolic(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);

void print_spec_interp_centroid(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);

void print_spec_interp_grandke(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);

void print_spec_interp_derivative(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);

void print_spec_interp_reassignment(Buffer* buf, const unsigned int spec_index,
																		unsigned int offset);
																		
void print_spec_interp_trigonometric(Buffer* buf, const unsigned int spec_index,
															unsigned int offset);

void print_spec_interp_charpentier(Buffer* buf, const unsigned int spec_index,
															unsigned int offset);

void print_spec_peaks(Buffer *buf, const unsigned int spec_index, unsigned int offset);

void print_spec_peaks_interp(Buffer *buf, const unsigned int spec_index,
				void (*interp_method)(Buffer* buf, const unsigned char channel,
					const unsigned int spec_index, const unsigned int bin_to_interp));

void print_spec_locmax(Buffer* buf, const unsigned int spec_index, unsigned int offset);

double trigonometric_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp);

double charpentier_interp(Buffer* buf, const unsigned char channel,
						const unsigned int spec_index, const unsigned int bin_to_interp);

void init_normalization(const Buffer *buf);


/*
void print_spec_interp(Buffer* buf, const unsigned int spec_index, unsigned int offset);
*/

#endif /* _SPEC_FUN_H_ */
