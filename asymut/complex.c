#include <math.h>
#include <stdio.h>
#include <string.h>

#include "trigonometric_constants.h"
#include "complex.h"

/********************************************************************************************/
/* C_CONJUGATE                                                                              */
/********************************************************************************************/
Complex c_conjugate(const Complex original)
{
	Complex conjugate;
	
	conjugate.re = original.re;
	conjugate.im = -original.im;
	
	return conjugate;
}

/********************************************************************************************/
/* C_POLAR_TO_COMPLEX                                                                       */
/********************************************************************************************/
Complex c_polar_to_complex(const double magnitude, const double phase)
{
	Complex c;
	
	c.re = magnitude * cos(phase);
	c.im = magnitude * sin(phase);
	
	return c;
}

/********************************************************************************************/
/* C_COMPLEX_TO_POLAR                                                                       */
/********************************************************************************************/
void c_complex_to_polar(const Complex c, double *mag, double *ph)
{
	*mag = c_modulus(c);
	*ph = c_phase(c);
}

/********************************************************************************************/
/* C_MODULUS_SQUARE                                                                         */
/********************************************************************************************/
double c_modulus_square(const Complex c)
{
    return c.re*c.re + c.im*c.im;
}

/********************************************************************************************/
/* C_MODULUS                                                                                */
/********************************************************************************************/
double c_modulus(const Complex c)
{
    return sqrt(c.re*c.re + c.im*c.im);
}

/********************************************************************************************/
/* C_SAFE_MODULUS                                                                           */
/********************************************************************************************/
/* supposedly numerically safer than c_modulus */
double c_safe_modulus(const Complex c)
{
	if (fabs(c.re) >= fabs(c.im)) {
		return fabs(c.re)*sqrt(1+(c.im*c.im)/(c.re*c.re));
	} else {
		return fabs(c.im)*sqrt(1+(c.re*c.re)/(c.im*c.im));	
	}
}

/********************************************************************************************/
/* C_PHASE                                                                                  */
/********************************************************************************************/
double c_phase(const Complex c)
{
    /* YES, it seems to work properly for the 2nd and 3rd quadrants */
    return atan2(c.im, c.re);
}

/********************************************************************************************/
/* C_COMPARE                                                                                */
/********************************************************************************************/
signed char c_compare(const Complex a, const Complex b)
{
    double x = c_modulus(a) - c_modulus(b);
    
    if (x)
        return (x > 0 ? 1 : -1);
        
    return 0;
}

/********************************************************************************************/
/* C_ARRAYS_COMPARE                                                                         */
/********************************************************************************************/
signed char c_arrays_compare(const Complex v1[], const Complex v2[], const unsigned long n)
{
    unsigned long i;

    for (i = 0; i < n; i++)
        if ( c_compare(v1[i], v2[i]) )
        	return ( c_compare(v1[i], v2[i]) );

    return 0;
}

/********************************************************************************************/
/* C_ADD                                                                                    */
/********************************************************************************************/
Complex c_add(const Complex a, const Complex b)
{
    Complex c;
    
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    
    return c;
}

/********************************************************************************************/
/* C_SUBTRACT                                                                               */
/********************************************************************************************/
Complex c_subtract(const Complex a, const Complex b)
{
    Complex c;
    
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    
    return c;

}    

/********************************************************************************************/
/* C_MULT                                                                                   */
/********************************************************************************************/
Complex c_mult(const Complex a, const Complex b)
{
    Complex c;
    
    c.re = (a.re * b.re) - (a.im * b.im);
    c.im = (a.re * b.im) + (a.im * b.re);
    
    return c;
}

/********************************************************************************************/
/* C_DIV                                                                                    */
/********************************************************************************************/
Complex c_div(const Complex a, const Complex b)
{
    Complex c;
    const double squared_b_modulus = (b.re*b.re + b.im*b.im);
    
    c.re = ((a.re * b.re) + (a.im * b.im))/squared_b_modulus;
    c.im = ((a.im * b.re) - (a.re * b.im))/squared_b_modulus;
    
    return c;
}


/********************************************************************************************/
/* C_DIV_SAFE                                                                               */
/********************************************************************************************/
Complex c_div_safe(const Complex a, const Complex b)
{
    Complex c;
    
    if (b.re >= b.im)
    {
    	double aux = b.im/b.re,
    			denominator = b.re + b.im*aux;

    	c.re = (a.re + a.im*aux)/denominator;
    	c.im = (a.im - a.re*aux)/denominator;
    }
    else
    {
    	double aux = b.re/b.im,
    			denominator = b.re*aux + b.im;
    			
		c.re = (a.re*aux + a.im)/denominator;
		c.im = (a.im*aux - a.re)/denominator;
    }
    
    return c;
}


/********************************************************************************************/
/* C_EXP                                                                                    */
/********************************************************************************************/
Complex c_exp(const double e)
{
    Complex c;
    
    c.re = cos(e);
    c.im = sin(e);
    
    return c;
}

/********************************************************************************************/
/* C_COPY                                                                                   */
/********************************************************************************************/
void c_copy(const Complex a, Complex* b)
{
    memcpy(b, &a, sizeof(Complex));
}

/********************************************************************************************/
/* REFERENCE_C_COPY                                                                         */
/********************************************************************************************/
void reference_c_copy(const Complex a, Complex* b)
{
    (*b).re = a.re;
    (*b).im = a.im;
}

/********************************************************************************************/
/* C_REAL_TO_COMPLEX                                                                        */
/********************************************************************************************/
Complex c_real_to_complex(const double real)
{
    Complex c;
    
    c.re = real;
    c.im = 0;
    
    return c;
}

/********************************************************************************************/
/* C_IMAGINARY_TO_COMPLEX                                                                   */
/********************************************************************************************/
Complex c_imaginary_to_complex(const double imaginary)
{
    Complex c;
    
    c.re = 0;
    c.im = imaginary;
    
    return c;
}

/********************************************************************************************/
/* C_DIV_BY_REAL                                                                            */
/********************************************************************************************/
Complex c_div_by_real(const Complex c, const double d)
{
    Complex y;
    
    y.re = c.re / d;
    y.im = c.im / d;
    
    return y;
}

/********************************************************************************************/
/* C_MULT_BY_REAL                                                                           */
/********************************************************************************************/
Complex c_mult_by_real(const Complex c, const double d)
{
    Complex y;

    y.re = c.re * d;
    y.im = c.im * d;

    return y;
}

/********************************************************************************************/
/* C_DISTANCE                                                                               */
/********************************************************************************************/
double c_distance(const Complex a, const Complex b)
{
    const double dr = (a.re - b.re), di = (a.im - b.im);
    
    return sqrt(dr*dr + di*di);
}

/********************************************************************************************/
/* C_PRINT                                                                                  */
/********************************************************************************************/
void c_print(const Complex c)
{
    printf("(%g, %g)\n", c.re, c.im);
}

/********************************************************************************************/
/* C_PRINT_THRESHOLDED                                                                      */
/********************************************************************************************/
void c_print_thresholded(const Complex c)
{
    c_print(c_threshold(c));
}

/********************************************************************************************/
/* C_ARRAY_PRINT                                                                            */
/********************************************************************************************/
void c_array_print(const Complex v[], const unsigned long n)
{
    unsigned long i;
    
    for (i = 0; i < n; i++)
        c_print(v[i]);
}

/********************************************************************************************/
/* C_ARRAY_PRINT_THRESHOLDED                                                                */
/********************************************************************************************/
void c_array_print_thresholded(const Complex v[], const unsigned long n)
{
    unsigned long i;

    for (i = 0; i < n; i++)
        c_print_thresholded(v[i]);
}

/********************************************************************************************/
/* C_ARRAY_MAG_PHASE_PRINT                                                                  */
/********************************************************************************************/
void c_array_mag_phase_print(const Complex v[], const unsigned long n)
{
    unsigned long i;
    
    for (i = 0; i < n; i++)
        printf("Mag: %lf, Pha: %lf * pi\n", c_modulus(v[i]), c_phase(v[i]) / PI);
}

/********************************************************************************************/
/* C_ARRAY_MAG_PHASE_PRINT_THRESHOLDED                                                      */
/********************************************************************************************/
void c_array_mag_phase_print_thresholded(const Complex v[], const unsigned long n)
{
    unsigned long i;

    for (i = 0; i < n; i++)
        printf("Mag: %lf, Pha: %lf * pi\n", c_modulus(c_threshold(v[i])),
        									c_phase(c_threshold(v[i])) / PI);
}

/********************************************************************************************/
/* C_THRESHOLD                                                                              */
/********************************************************************************************/
Complex c_threshold(const Complex c)
{
    Complex thresholded;

    thresholded.re = (fabs(c.re) < THRESHOLD_VALUE ? 0 : c.re);
    thresholded.im = (fabs(c.im) < THRESHOLD_VALUE ? 0 : c.im);

    return thresholded;
}
