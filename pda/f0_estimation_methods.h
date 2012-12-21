#ifndef _F0_ESTIMATION_METHODS_H_
#define _F0_ESTIMATION_METHODS_H_

#include "data_structures.h"

double min(double *v, const unsigned int n);

double max(double *v, const unsigned int n);

int round_to_int(const double x);

unsigned int max_index(double *v, const unsigned int n);

unsigned int threshold(double *v, const unsigned int n, double threshold_value);

unsigned int threshold2(double *v, double *w, const unsigned int n, double threshold_value);

double sum(double *v, const unsigned int n);

double aud_sum(double *v, const unsigned int n);

double average(double *v, const unsigned int n);

double average2(double *v, const unsigned int n);

double odd_average2(double *v, const unsigned int n);

double variance(double *v, const unsigned int n);

double variance2(double *v, const unsigned int n);

double weighted_variance(double *v, double *w, const unsigned int n);

double lin_reg_f0_est(double freq[], double harm[], const unsigned int n);

double average_f0_est(double freq[], double amp[], const unsigned int n);

double cb_average_f0_est(double freq[], double amp[], const unsigned int n,
																	const double cb_weight);

double f0_error(double f0, double freq[], double amp[], const unsigned int n);

double f0_cb_error(double f0, double freq[], double amp[], const unsigned int n,
																	const double cb_weight);

void spec_smooth(double freq[], double amp[], const unsigned int n);

double full_average_f0_est(double freq[], double amp[], const unsigned int n,
	const double min_f0, const double max_f0,
    	const double cb_weight_error, const double cb_weight_f0_est);
    	
double smoothed_average_f0_est(double freq[], double amp[], const unsigned int n,
	const double min_f0, const double max_f0,
    	const double cb_weight_error, const double cb_weight_f0_est);

double mitre_f0_est(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks,
	const double time, const EstimationParameters *estPar);

double erb_f0_est(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks,
	const double time, const EstimationParameters *estPar);

double getFreqByAmpQuantile(double peaksFreq[], double peaksAmp[], const unsigned int numOfPeaks, const double quantile);


double absoluteRelativeError(const double f1, const double f2);

double averageAbsF0Error(double f0, double freq[], double amp[], const double total_energy, const unsigned int n);

void prepare_critical_band_sum(unsigned int n);

void finalize_critical_band_sum(void);

double critical_band_sum(double *v, const unsigned int n);

double new_critical_band_sum(const double f0, double freq[], double mag[], const unsigned int n);

double erb_sum(double freq[], double mag[], const unsigned int n);

double bark_sum(double freq[], double mag[], const unsigned int n);

double critical_band_average2(double *v, const unsigned int n);

double odd_critical_band_average2(double *v, const unsigned int n);

#endif
