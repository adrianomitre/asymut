#include "data_structures.h"

signed char parse_arguments(const int n_args, char** argument,
							Audio_Stream* input, FFT_properties* fft);

void print_current_config(Buffer *buf);

unsigned char is_power_of_2(unsigned long n);
unsigned char verify_extension (const char* name, const char* extension);
signed char open_audiofile(Audio_Stream* stream);
signed char close_audiofile(Audio_Stream* stream);
short read_wav_header(Audio_Stream* stream);

void print_fft_prop_state(FFT_properties* fft_p, const int n);

float c8_to_float(const unsigned char sample);
float c16_to_float(const unsigned char sample[2]);

unsigned short round_double_to_ushort(double val);

double princ_det(double angle);

void* emalloc(const size_t size);
void* ecalloc(const size_t size);

void calc_mov_median(double signal[], unsigned int sig_len, unsigned short win_len,
							unsigned short med_index, double scaling_factor,
       						double const_part, double median[]);

void selection_sort(double v[], unsigned int n);

double log2(double x);

double hertz_to_bark(const double f);
double erb_to_hertz(const double b);
double hertz_to_bark(const double f);

#define max(a, b) b > a ? b : a
#define min(a, b) b < a ? b : a
