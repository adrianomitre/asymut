#ifndef _AUX_FUN_H_
#define _AUX_FUN_H_

signed char parse_arguments(const int n_args, char** argument, FILE **F0_fp, FILE **onset_fp);
void* emalloc(const size_t size);
void* ecalloc(const size_t size);
double log2(double x);

void* erealloc(void *p, const size_t size);
void free_2d(void **p, const unsigned int n);

double hertz_to_bark(const double f);
double hertz_to_erb(const double f);


#endif /* _AUX_FUN_H_ */
