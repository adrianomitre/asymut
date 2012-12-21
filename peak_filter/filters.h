#ifndef _FILTERS_H_
#define _FILTERS_H_

double min(double *v, const unsigned int n);

double max(double *v, const unsigned int n);

int round_to_int(const double x);

unsigned int get_max_index(double *v, const unsigned int n);

unsigned int threshold(double *v, const unsigned int n, double threshold_value);

unsigned int threshold2(double *v, double *w, const unsigned int n, double threshold_value);

double sum(double *v, const unsigned int n);

double average(double *v, const unsigned int n);

double cubic_average(double *v, const unsigned int n);

double average2(double *v, const unsigned int n);

double variance(double *v, const unsigned int n);

double variance2(double *v, const unsigned int n);

double weighted_variance(double *v, double *w, const unsigned int n);

unsigned int remove_masked(double amp[], double freq[], const unsigned int n,
											const double low_ratio, const double high_ratio);

unsigned int remove_spurious_peaks(double amp[], double freq[], const unsigned int n);

unsigned int new_filter(double *amp, double *freq, const unsigned int n,
																double max_per_octave_decay);

unsigned int amg_filter(double *amp, double *freq, const unsigned int n,
																double max_adj_harm_gap);

unsigned int multi_filter(double *mag, double *freq, const unsigned int num_of_peaks,
							const double min_freq, const double max_freq,
							const double min_mag,
							const double max_mag_gap,
							const double max_adj_mag_gap);

#endif /* _FILTERS_H_ */
