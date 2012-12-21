#ifndef _COMPLEX_H_
#define _COMPLEX_H_

static const double THRESHOLD_VALUE = 0.000000000001;

typedef struct complex_t {
    double re;
    double im;
} Complex;

Complex c_conjugate(const Complex original);

double c_modulus(const Complex c);
double c_safe_modulus(const Complex c);
double c_modulus_square(const Complex c);
double c_phase(const Complex c);

signed char c_compare(const Complex a, const Complex b);
signed char c_arrays_compare(const Complex v1[], const Complex v2[], const unsigned long n);

Complex c_add(const Complex a, const Complex b);
Complex c_subtract(const Complex a, const Complex b);
Complex c_mult(const Complex a, const Complex b);
Complex c_div(const Complex a, const Complex b);
Complex c_exp(const double e);
Complex c_conjug(const Complex c);
Complex c_polar_to_complex(const double magnitude, const double phase);
void c_complex_to_polar(const Complex c, double *mag, double *ph);

void c_copy(const Complex a, Complex* b);

Complex c_real_to_complex(const double r);
Complex c_imaginary_to_complex(const double i);
Complex c_div_by_real(const Complex c, const double d);
Complex c_mult_by_real(const Complex c, const double d);

double c_distance(const Complex a, const Complex b);

void c_print(const Complex c);
void c_print_thresholded(const Complex c);
void c_array_print(const Complex v[], const unsigned long n);
void c_array_print_thresholded(const Complex v[], const unsigned long n);
void c_array_mag_phase_print(const Complex v[], const unsigned long n);
void c_array_mag_phase_print_thresholded(const Complex v[], const unsigned long n);

Complex c_threshold(const Complex c);

/* no longer being used, for the sake of computational speed */
void reference_c_copy(const Complex a, Complex* b);

#endif
