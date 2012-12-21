#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "numeric_constants.h"
#include "trigonometric_constants.h"
#include "config_constants.h"
#include "config_variables.h"
#include "complex.h"
#include "aux_fun.h"
#include "stdft.h"

static unsigned int *rev_table;
static double *sin_table, *cos_table, *unpack_sin_table, *unpack_cos_table, *out_re, *out_im;
static Complex *packed, *unpacked;

/********************************************************************************************/
/* REFERENCE_DFT                                                                            */
/********************************************************************************************/
void reference_dft(const Complex s[], Complex ft[], const unsigned long n)
{
	if (n%2) {
		fprintf(stderr, "Error: DFT is not defined for odd-sized window.\n");
		return;
	}

	const long half_n = n/2;
    long t, w;
    unsigned long t_ind, w_ind;

    for (w = -half_n+1; w <= half_n; w++) {
    	w_ind = w+half_n-1;
    	ft[w_ind].re = ft[w_ind].im = 0;
	    for (t = -half_n; t < half_n; t++) {
	    	t_ind = t + half_n;
	    	c_copy(c_add(c_mult(s[t_ind], c_exp(-TWO_PI*t*w/(double)n)), ft[w_ind]), &ft[w_ind]);
	    }
    }
}

/********************************************************************************************/
/* UNSIGNED_INDEX                                                                           */
/********************************************************************************************/
unsigned long unsigned_index(const signed long n, const unsigned long window_size)
{
    unsigned long i = n + (window_size / 2);

    if (i >= 0) return i;
    else
    {
        fprintf(stderr, "Serious index error! Returning first valid array index "
        				"(i.e., zero)...\n");
        return 0;
    }
}

/********************************************************************************************/
/* PREPARE_FFT_REAL_SIGNAL_TO_MAG_PHASE                                                     */
/********************************************************************************************/
void old_prepare_fft_real_signal_to_mag_phase(FFT_properties* fft_properties)
{
  const unsigned int n = fft_properties->window.length * fft_properties->zero_padding_ratio;

  out_re = (double*)emalloc(n * sizeof(double));
  out_im = (double*)emalloc(n * sizeof(double));

  init_rev_table(n);
  init_sin_table(n);
  init_cos_table(n);
}

/********************************************************************************************/
/* FINALIZE_FFT_REAL_SIGNAL_TO_MAG_PHASE                                                    */
/********************************************************************************************/
void old_finalize_fft_real_signal_to_mag_phase(void)
{
  free_rev_table();
  free_sin_table();
  free_cos_table();

  if (out_re) free(out_re);
  if (out_im) free(out_im);
}

/********************************************************************************************/
/* OLD_FFT_REAL_SIGNAL_TO_MAG_PHASE                                                         */
/********************************************************************************************/
double old_fft_real_signal_to_mag_phase(double signal[], FFT_properties* fft_properties,
														double magnitude[], double phase[])
{
    const unsigned int core_len = fft_properties->window.length,
    					n = core_len * fft_properties->zero_padding_ratio;

    double w_m_re, w_re, t_re, u_re,
    	   w_m_im, w_im, t_im, u_im;

    unsigned int k, m, j; /* try using the 'register' qualifier */

    for (k = 0; k < core_len; k++)
    {
        out_re[rev_table[k]] = fft_properties->window.content[k] * signal[k];
        out_im[k] = 0;
    }
    for (k = core_len; k < n; k++)
    {
        out_re[rev_table[k]] = 0;
        out_im[k] = 0;
    }


    for (m = 2; m <= n; m = (m << 1))
    {
    	w_m_re = cos_table[m];
     	w_m_im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w_re = 1;
            w_im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
            	t_re = (w_re * out_re[k + j + (m >> 1)]) - (w_im * out_im[k + j + (m >> 1)]);
             	t_im = (w_re * out_im[k + j + (m >> 1)]) + (w_im * out_re[k + j + (m >> 1)]);

             	u_re = out_re[k + j];
              	u_im = out_im[k + j];

              	out_re[k + j] = u_re + t_re;
              	out_im[k + j] = u_im + t_im;

              	out_re[k + j + (m >> 1)] = u_re - t_re;
              	out_im[k + j + (m >> 1)] = u_im - t_im;

                t_re = (w_re * w_m_re) - (w_im * w_m_im);
                t_im = (w_re * w_m_im) + (w_im * w_m_re);
                w_re = t_re;
                w_im = t_im;
            }
        }
    }

    for (k = 0, m = n >> 1; k < m; k++)
    {
        t_re = out_re[k + 1];
        t_im = out_im[k + 1];

        magnitude[k] = sqrt(t_re*t_re + t_im*t_im);
        phase[k] = atan2(t_im, t_re);
    }
    magnitude[m - 1] /= 2;

    t_re = out_re[0];
    t_im = out_im[0];

    return sqrt(t_re*t_re + t_im*t_im);
}

/********************************************************************************************/
/* FFT_REAL_SIGNAL_TO_MAG                                                                   */
/********************************************************************************************/
void fft_real_signal_to_mag(double signal[], FFT_properties* fft_properties,
																		double magnitude[])
{
    static unsigned char first_call = TRUE;
    static double *out_re, *out_im;

    double w_m_re, w_re, t_re, u_re,
    	   w_m_im, w_im, t_im, u_im;

    unsigned int k, m, j; /* try using the 'register' qualifier */

    if (fft_properties == NULL)
    {
        if (!first_call)
        {
        	if (out_re) free(out_re);
        	if (out_im) free(out_im);

        	first_call = TRUE;
        }

        return;
    }


   {
    const unsigned int core_len = fft_properties->window.length,
    					n = fft_properties->window.length * fft_properties->zero_padding_ratio;

    if (first_call == TRUE)
    {
        out_re = (double*)emalloc(n * sizeof(double));
        out_im = (double*)emalloc(n * sizeof(double));

        first_call = FALSE;
    }

    for (k = 0; k < core_len; k++)
    {
        out_re[rev_table[k]] = fft_properties->window.content[k] * signal[k];
        out_im[k] = 0;
    }
    for (k = core_len; k < n; k++)
    {
        out_re[k] = 0;
        out_im[k] = 0;
    }


    for (m = 2; m <= n; m = (m << 1))
    {
    	w_m_re = cos_table[m];
     	w_m_im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w_re = 1;
            w_im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
            	t_re = (w_re * out_re[k + j + (m >> 1)]) - (w_im * out_im[k + j + (m >> 1)]);
             	t_im = (w_re * out_im[k + j + (m >> 1)]) + (w_im * out_re[k + j + (m >> 1)]);

             	u_re = out_re[k + j];
              	u_im = out_im[k + j];

              	out_re[k + j] = u_re + t_re;
              	out_im[k + j] = u_im + t_im;

              	out_re[k + j + (m >> 1)] = u_re - t_re;
              	out_im[k + j + (m >> 1)] = u_im - t_im;

                t_re = (w_re * w_m_re) - (w_im * w_m_im);
                t_im = (w_re * w_m_im) + (w_im * w_m_re);
                w_re = t_re;
                w_im = t_im;
            }
        }
    }

    for (k = 0, m = n >> 1; k < m; k++)
    {
        t_re = out_re[k + 1];
        t_im = out_im[k + 1];

        magnitude[k] = sqrt(t_re*t_re + t_im*t_im);
    }
    magnitude[m - 1] /= 2;
   }
}

/********************************************************************************************/
/* PURELY_REAL_SIGNAL_FFT                                                                   */
/********************************************************************************************/
/* efficiently computes the FFT of a purely real signal to a previously allocated pair
 * of double-type arrays */
void purely_real_signal_fft(const double signal[], FFT_properties* fft_properties,
																	double mag[], double ph[])
{
	const unsigned int effective_fft_size = fft_properties->window.length/2,
				total_fft_size = effective_fft_size * fft_properties->zero_padding_ratio;
	unsigned int i;
	Complex *in, *out, spec, spec_even, spec_odd, minus_j;

	in = (Complex*)emalloc(total_fft_size*sizeof(Complex));
	out = (Complex*)emalloc(total_fft_size*sizeof(Complex));

	minus_j.re = 0;
	minus_j.im = -1;

	for (i = 0; i < effective_fft_size; i++) {
		in[i].re = signal[2*i] * fft_properties->window.content[2*i];
		in[i].im = signal[2*i + 1] * fft_properties->window.content[2*i + 1];
	}
	for (i = effective_fft_size; i < total_fft_size; i++) {
		in[i].re = 0; /* this must not be done here, every call, for speed's sake */
		in[i].im = 0; /* but rather on the same place where in[] will be allocated */
	}

	fft(in, out, total_fft_size);

	for (i = 0; i < total_fft_size; i++) {
		c_copy(c_div_by_real(c_add(out[i], c_conjugate(out[total_fft_size-i])), 2), &spec_even);
		c_copy(c_mult(c_div_by_real(c_subtract(out[i], c_conjugate(out[total_fft_size-i])), 2), minus_j), &spec_odd);
		c_copy(c_add(spec_even, c_mult(spec_odd, c_exp(-TWO_PI*i/total_fft_size))), &spec);

		if (mag) mag[i] = c_modulus(spec); /* this IF slow down perceively? */
		if (ph) ph[i] = c_phase(spec); /* this IF slow down perceively? */
	}
}

/********************************************************************************************/
/* FFT                                                                                      */
/********************************************************************************************/
void fft(const Complex in[], Complex out[], const unsigned int n)
{
    Complex w_m, w, t, u; /* try using the 'register' qualifier */
    unsigned int k, m, j; /* try using the 'register' qualifier */

    for (k = 0; k < n; k++)
    {
        out[rev_table[k]].re = in[k].re;
        out[rev_table[k]].im = in[k].im;
    }

    for (m = 2; m <= n; m = (m << 1))
    {
    	w_m.re = cos_table[m];
     	w_m.im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w.re = 1;
            w.im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
            	t.re = (w.re * out[k + j + (m >> 1)].re) - (w.im * out[k + j + (m >> 1)].im);
             	t.im = (w.re * out[k + j + (m >> 1)].im) + (w.im * out[k + j + (m >> 1)].re);

             	u.re = out[k + j].re;
              	u.im = out[k + j].im;

              	out[k + j].re = u.re + t.re;
              	out[k + j].im = u.im + t.im;

              	out[k + j + (m >> 1)].re = u.re - t.re;
              	out[k + j + (m >> 1)].im = u.im - t.im;

                t.re = (w.re * w_m.re) - (w.im * w_m.im);
                t.im = (w.re * w_m.im) + (w.im * w_m.re);
                w.re = t.re;
                w.im = t.im;
            }
        }
    }
}


/********************************************************************************************/
/* PREVIOUS_FFT                                                                             */
/********************************************************************************************/
void previous_fft(const Complex in[], Complex out[], const unsigned int n)
{
    Complex w_m, w, t, u; /* try using the 'register' qualifier */
    unsigned int k, m, j; /* try using the 'register' qualifier */

    for (k = 0; k < n; k++)
    {
        out[rev_table[k]].re = in[k].re; /* TODO: what about rev_table_k++? */
        out[k].im = 0; /* this is possible because the signal is entirely real*/
    }

    for (m = 2; m <= n; m = (m << 1))
    {
    	w_m.re = cos_table[m];
     	w_m.im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w.re = 1;
            w.im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
            	t.re = (w.re * out[k + j + (m >> 1)].re) - (w.im * out[k + j + (m >> 1)].im);
             	t.im = (w.re * out[k + j + (m >> 1)].im) + (w.im * out[k + j + (m >> 1)].re);

             	u.re = out[k + j].re;
              	u.im = out[k + j].im;

              	out[k + j].re = u.re + t.re;
              	out[k + j].im = u.im + t.im;

              	out[k + j + (m >> 1)].re = u.re - t.re;
              	out[k + j + (m >> 1)].im = u.im - t.im;

                t.re = (w.re * w_m.re) - (w.im * w_m.im);
                t.im = (w.re * w_m.im) + (w.im * w_m.re);
                w.re = t.re;
                w.im = t.im;
            }
        }
    }
}


/********************************************************************************************/
/* FFT2                                                                                     */
/********************************************************************************************/
void fft2(const Complex in[], Complex out[], const unsigned int n)
{
    Complex w_m, w, t, u; /* try using the 'register' qualifier */
    unsigned int k, m, j; /* try using the 'register' qualifier */

    for (k = 0; k < n; k++)
    {
        out[rev_table[k]].re = in[k].re; /* TODO: what about rev_table_k++? */
        out[k].im = 0; /* this is possible because the signal is entirely real*/
    }

    for (m = 2; m <= n; m = (m << 1))
    {
    	w_m.re = cos_table[m];
     	w_m.im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w.re = 1;
            w.im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
            	t.re = (w.re * out[k + j + (m >> 1)].re) - (w.im * out[k + j + (m >> 1)].im);
             	t.im = (w.re * out[k + j + (m >> 1)].im) + (w.im * out[k + j + (m >> 1)].re);

             	u.re = out[k + j].re;
              	u.im = out[k + j].im;

              	out[k + j].re = u.re + t.re;
              	out[k + j].im = u.im + t.im;

              	out[k + j + (m >> 1)].re = u.re - t.re;
              	out[k + j + (m >> 1)].im = u.im - t.im;

                t.re = (w.re * w_m.re) - (w.im * w_m.im);
                t.im = (w.re * w_m.im) + (w.im * w_m.re);
                w.re = t.re;
                w.im = t.im;
            }
        }
    }
}

/********************************************************************************************/
/* REFERENCE_FFT                                                                            */
/********************************************************************************************/
void reference_fft(const Complex in[], Complex out[], const unsigned int n)
{
    Complex w_m, w, t, u;
    unsigned int k, m, j;

    for (k = 0; k < n; k++)
        c_copy(in[k], &out[rev(k, n)]);

    for (m = 2; m <= n; m = (m << 1))
    {
        c_copy(c_exp(TWO_PI / m), &w_m);
        for (k = 0; k < n; k += m)
        {
            c_copy(c_real_to_complex(1.0), &w);
            for (j = 0; j < (m >> 1); j++)
            {
                c_copy(c_mult(w, out[k + j + (m >> 1)]), &t);
                c_copy(out[k + j], &u);
                c_copy(c_add(u, t), &out[k + j]);
                c_copy(c_subtract(u, t), &out[k + j + (m >> 1)]);
                c_copy(c_mult(w, w_m), &w);
            }
        }
    }

    /* TODO: reimplement without using auxiliary array */
	Complex *aux = (Complex*)emalloc(n*sizeof(Complex));
	memmove(aux, out, n*sizeof(Complex));
	for (k = 0; k < n; k++) c_copy(aux[k], &out[(k+n/2-1)%n]);
	free(aux);
}

/********************************************************************************************/
/* REFERENCE_FFT2                                                                           */
/********************************************************************************************/
void reference_fft2(const Complex in[], Complex out[], const unsigned short n)
{
    Complex w_m, w, t, u;
    unsigned short k, m, j;

    for (k = 0; k < n; k++)
        c_copy(in[k], &out[rev(k, n)]);

    for (m = 2; m <= n; m = (m << 1))
    {
        c_copy(c_exp(TWO_PI / m), &w_m);
        for (k = 0; k < n; k += m)
        {
            c_copy(c_real_to_complex(1.0), &w);
            for (j = 0; j < (m >> 1); j++)
            {
                c_copy(c_mult(w, out[k + j + (m >> 1)]), &t);
                c_copy(out[k + j], &u);
                c_copy(c_add(u, t), &out[k + j]);
                c_copy(c_subtract(u, t), &out[k + j + (m >> 1)]);
                c_copy(c_mult(w, w_m), &w);
            }
        }
    }
}


/********************************************************************************************/
/* INIT_REV_TABLE                                                                           */
/********************************************************************************************/
void init_rev_table(const unsigned int n)
{
    unsigned int k, i, j, half_n = n >> 1, reverse;

    rev_table = (unsigned int*)emalloc(n * sizeof(unsigned int));

    for (k = 0; k < n; k++)
    {
        for (i = 1, j = half_n, reverse = 0; j > 0; i = (i << 1), j = (j >> 1))
            if (k & j)
                reverse += i;
        rev_table[k] = reverse;
    }
}

/********************************************************************************************/
/* FREE_REV_TABLE                                                                           */
/********************************************************************************************/
void free_rev_table(void)
{
    if (rev_table) free(rev_table);
}

/********************************************************************************************/
/* INIT_SIN_TABLE                                                                           */
/********************************************************************************************/
void init_sin_table(const unsigned int n)
{
    unsigned int m;

    sin_table = (double*)emalloc((n + 1) * sizeof(double));

	/*
		for (m = 0; m <= n; m++)
			sin_table[m] = 0;
	*/

    for (m = 2; m <= n; m <<= 1)
        sin_table[m] = sin(TWO_PI / m);

	/*
		const unsigned int half_n = n/2;
		for (m = 1; m <= half_n; m++)
		{
			double aux = sin_table[m];
			sin_table[m] = sin_table[m + half_n];
			sin_table[m + half_n] = aux;
		}
	*/
}

/********************************************************************************************/
/* FREE_SIN_TABLE                                                                           */
/********************************************************************************************/
void free_sin_table(void)
{
    if (sin_table) free(sin_table);
}

/********************************************************************************************/
/* INIT_COS_TABLE                                                                           */
/********************************************************************************************/
void init_cos_table(const unsigned int n)
{
    unsigned int m;

    cos_table = (double*)emalloc((n + 1)* sizeof(double));

	/*
		for (m = 0; m <= n; m++)
			cos_table[m] = 0;
	*/

    for (m = 2; m <= n; m <<= 1)
        cos_table[m] = cos(TWO_PI / m);

	/*
		const unsigned int half_n = n/2;
		for (m = 1; m <= half_n; m++)
		{
			double aux = sin_table[m];
			cos_table[m] = cos_table[m + half_n];
			cos_table[m + half_n] = aux;
		}
	*/
}

/********************************************************************************************/
/* FREE_COS_TABLE                                                                           */
/********************************************************************************************/
void free_cos_table(void)
{
    if (cos_table) free(cos_table);
}

/********************************************************************************************/
/* REV                                                                                      */
/********************************************************************************************/
unsigned int rev(unsigned int k, unsigned int n)
{
    unsigned int i, j, half_n = n >> 1, reverse = 0;

    for (i = 1, j = half_n; j > 0; i = (i << 1), j = (j >> 1))
        if (k & j)
            reverse += i;

    return reverse;
}

/********************************************************************************************/
/* INIT_WINDOW                                                                              */
/********************************************************************************************/
void init_window(Window *w)
{
	void (*APOD_FUN_PTR)(Window *w);

    if (!is_power_of_2(w->length))
    {
   	    fprintf(stderr, "Error in init_window(): window length size must be a power of 2.\n");
   	    exit ( ERROR );
    }

	w->content = (double*)emalloc(w->length * sizeof(double));
	w->derivative = (double*)emalloc(w->length * sizeof(double));

    switch (w->apodization_function.type)
    {
        case TRIANGULAR:
            APOD_FUN_PTR = triangular;
            break;
        case RECTANGULAR:
            APOD_FUN_PTR = rectangular;
            break;
        case BLACKMAN_HARRIS:
            APOD_FUN_PTR = blackman_harris;
            break;
        case FLAT_TOP:
            APOD_FUN_PTR = flat_top;
            break;
        case HANN:
            APOD_FUN_PTR = hann;
            break;
        case BLACKMAN:
            APOD_FUN_PTR = blackman;
            break;
        case NUTTALL3:
            APOD_FUN_PTR = nuttall3;
            break;
        case NUTTALL4:
            APOD_FUN_PTR = nuttall4;
            break;
        case NUTTALL5:
            APOD_FUN_PTR = nuttall5;
            break;
        case NUTTALL6:
            APOD_FUN_PTR = nuttall6;
            break;
        case NUTTALL7:
            APOD_FUN_PTR = nuttall7;
            break;
        case NUTTALL8:
            APOD_FUN_PTR = nuttall8;
            break;
        case NUTTALL9:
            APOD_FUN_PTR = nuttall9;
            break;
        case NUTTALL10:
            APOD_FUN_PTR = nuttall10;
            break;
        case NUTTALL11:
            APOD_FUN_PTR = nuttall11;
            break;
        case NUTTALL12:
            APOD_FUN_PTR = nuttall12;
            break;
        case HAMMING:
            APOD_FUN_PTR = hamming;
            break;
        case NUTTALL14:
            APOD_FUN_PTR = nuttall14;
            break;
        case NUTTALL15:
            APOD_FUN_PTR = nuttall15;
            break;
        case GAUSSIAN:
        	APOD_FUN_PTR = gaussian;
        	break;
        case HANNING_POISSON:
        	APOD_FUN_PTR = hanning_poisson;
        	break;
        case HELIE_A_W1:
        	APOD_FUN_PTR = helie_a_w1;
        	break;
        case HELIE_A_W6:
        	APOD_FUN_PTR = helie_a_w6;
        	break;
        case MOD_BARLETT_HANN:
        	APOD_FUN_PTR = mod_barlett_hann;
        	break;
        case BOHMAN:
        	APOD_FUN_PTR = bohman;
        	break;
        default:
        	fprintf(stderr, "Fatal error in init_window(): unexpected apodization function.\n");
        	exit(1);
    }

    APOD_FUN_PTR(w);
}

/********************************************************************************************/
/* PRINT_WINDOW_TIME_RESPONSE                                                               */
/********************************************************************************************/
void print_window_time_response(Window* w)
{
    const unsigned int length = w->length, half_length = length/2;
    unsigned int i;

    for (i = 0; i < length; i++)
    	printf("%g %g\n", ((double)i - half_length + 0.5) / half_length, w->content[i]);
}

/********************************************************************************************/
/* PRINT_WINDOW_FREQ_RESPONSE                                                               */
/********************************************************************************************/
void print_window_freq_response(FFT_properties *fft_p, unsigned int max_bin)
{
    unsigned int i;
    double *mag, *ph, *signal;

    if (!fft_p) return;

    {
    const unsigned short zero_padding_ratio = fft_p->zero_padding_ratio;
    const unsigned int core_len = fft_p->window.length,
    	    			total_len = core_len * zero_padding_ratio,
            			num_of_bins = total_len / 2;

    if (max_bin == 0 || max_bin > num_of_bins)
    	max_bin = num_of_bins / zero_padding_ratio;

    signal = (double*)emalloc(core_len * sizeof(double));
    mag = (double*)emalloc((num_of_bins + 1) * sizeof(double));
    ph = (double*)emalloc(num_of_bins * sizeof(double));

    for (i = 0; i < core_len; i++)
    	signal[i] = 1;


    old_prepare_fft_real_signal_to_mag_phase(fft_p);

    mag[0] = old_fft_real_signal_to_mag_phase(signal, fft_p, mag + 1, ph);

    old_finalize_fft_real_signal_to_mag_phase();

    for (i = max_bin * zero_padding_ratio; i > 0; i--)
    {
        if (mag[i] > 0) printf("%g %g\n",  (-1) * (double)i / zero_padding_ratio, 20*log10(mag[i]/mag[0]) );
    }
    for (i = 0; i <= max_bin * zero_padding_ratio; i++)
    {
        if (mag[i] > 0) printf("%g %g\n", (double)i / zero_padding_ratio, 20*log10(mag[i]/mag[0]) );
    }

    free(signal);
    free(mag);
    free(ph);
    }
}

/********************************************************************************************/
/* FREE_WINDOW                                                                              */
/********************************************************************************************/
void free_window(Window* w)
{
    if (w->content) free(w->content);
}

/********************************************************************************************/
/* TRIANGULAR                                                                               */
/********************************************************************************************/
void triangular(Window* w)
{
    const unsigned int length = w->length;
    const double scale = 2.0/length;
    double t;
    unsigned int i;

    w->total_area = 0;
    for (i = 0, t = length/-2.0; i < length; i++, t+=1)
    {
     	w->content[i] = 1.0 - fabs(t) * scale;
    	w->total_area += w->content[i];
    }
}

/********************************************************************************************/
/* OLD_TRIANGULAR                                                                           */
/********************************************************************************************/
/* same problem as old_raised_cosine */
void old_triangular(Window* w)
{
    const unsigned int length = w->length, half_length = length/2;
    const double scale = (double)1/half_length;
    unsigned int i;

    w->total_area = 0;
    for (i = 0; i < half_length; i++)
    {
     	w->content[i] = (i + 1) * scale;
    	w->total_area += w->content[i];
    }
    for (i = half_length; i < length; i++)
    {
    	w->content[i] = (length - i) * scale;
    	w->total_area += w->content[i];
    }
}


/********************************************************************************************/
/* RECTANGULAR                                                                              */
/********************************************************************************************/
void rectangular(Window* w)
{
    const unsigned int num_of_coef = 1;
    const double coef[1] = {1};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* BLACKMAN_HARRIS                                                                          */
/********************************************************************************************/
void blackman_harris(Window* w)
{
    const unsigned char num_of_coef = 4;
	const double coef[4] = {0.35875, 0.48829, 0.14128, 0.01168};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* FLAT_TOP                                                                                 */
/********************************************************************************************/
void flat_top(Window* w)
{
    const unsigned char num_of_coef = 5;
    /* 
     * References:
     * [1] ISO 18431-2:2004, Mechanical vibration and shock -- Signal
     *       processing -- Part 2: Time domain windows for Fourier Transform analysis
     * [2] Cortes et al, "A new class of flat-top windows for exposure assessment in magnetic
             field measurements", 2007
     * [3] Tran et al, "Window Design and Enhancement using Chebyshev Optimization", 2004
     *
     * According to [1], the 5-term flattop coefficients are
     *   {1, 1.933, 1.286, 0.388, 0.0322}
     * Among the flattops presented [2], the one I found most interesting is the
     *   fast-decaying / minimum-sidelobe 5-term, with 30 dB/oct asymptotic decay and maximal
     *   sidelobel level at -79.06 dB, whose coefficients are
     *   {0.20142488, 0.39291808, 0.28504554, 0.10708192, 0.01352957}. This was the chosen one.
     * Somewhere I found that the 3-term flattop coefficients are
     *   {0.2810639, 0.5208972, 0.1980399}
     */
    const double coef[5] = {0.20142488, 0.39291808, 0.28504554, 0.10708192, 0.01352957};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* FLAT_TOP_MATLAB                                                                          */
/********************************************************************************************/
/* MATLAB coefficients, as reported in [1], are are {0.21557895, 0.41663158, 0.277263158,
   0.083578947, 0.006947368}, but they sum to 1.000000003. The values used below were
   normalized to sum to unity.
   
   [1] http://www.mathworks.com/access/helpdesk/help/toolbox/signal/index.html?/access/helpdesk/help/toolbox/signal/flattopwin.html&http://www.google.com.br/search?q=flat+top+window&ie=utf-8&oe=utf-8&aq=t&rls=com.ubuntu:pt-BR:unofficial&client=firefox-a
 */
/* TODO: make it callable */
void flat_top_matlab(Window* w)
{
    const unsigned char num_of_coef = 5;
    const double coef[5] = {0.215578949353263, 0.416631578750105, 0.277263157168211,
                            0.0835789467492632, 0.0069473679791579};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* HANN                                                                                     */
/********************************************************************************************/
void hann(Window* w)
{
    const unsigned int num_of_coef = 2;
    const double coef[2] = {0.5, 0.5};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* BLACKMAN                                                                                 */
/********************************************************************************************/
void blackman(Window* w)
{
    const unsigned int num_of_coef = 3;
    const double coef[3] = {0.42, 0.5, 0.08};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL3  - a.k.a. "exact" Blackman                                                      */
/********************************************************************************************/
void nuttall3(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {(double)7938/18608, (double)9240/18608, (double)1430/18608};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL4                                                                                 */
/********************************************************************************************/
void nuttall4(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {0.42323, 0.49755, 0.07922};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL5                                                                                 */
/********************************************************************************************/
void nuttall5(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {0.44959, 0.49364, 0.05677};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL6                                                                                 */
/********************************************************************************************/
void nuttall6(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {0.35875, 0.48829, 0.14128, 0.01168};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL7                                                                                 */
/********************************************************************************************/
void nuttall7(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {0.40217, 0.49703, 0.09892, 0.00188};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL8                                                                                 */
/********************************************************************************************/
void nuttall8(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {0.375, 0.5, 0.125};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL9                                                                                 */
/********************************************************************************************/
void nuttall9(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {0.40897, 0.5, 0.09103};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL10                                                                                */
/********************************************************************************************/
void nuttall10(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {(double)10/32, (double)15/32, (double)6/32, (double)1/32};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL11                                                                                */
/********************************************************************************************/
void nuttall11(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {0.338946, 0.481973, 0.161054, 0.018027};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL12                                                                                */
/********************************************************************************************/
void nuttall12(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {0.355768, 0.487396, 0.144232, 0.012604};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* HAMMING                                                                                  */
/********************************************************************************************/
/* References:
 *  [1] Harris, "Use of Windows for Harmonic Analysis", 78
 *  [2] Nuttall, "Some Windows with Very Good Sidelobe Behavior", 81
 *  [3] Albrecht, "A family of cosine-sum windows for high-resolution measurements", 2001
 *  [4] Kulkarni and Lahiri, "Improved Sidelobe Performance of Cosine Series Functions", 99
 *
 * Peak sidelobe is -43 dB for alpha 0.54 [1], -43.19 for alpha 0.53836 [2],
 * -43.187 for 5.383553946707251e-001 [3], and -47.4 dB for alpha 0.54576 [4],
 * but not sure if it is the same window in [4]. Asymptotic decay is 6 dB/octave */
void hamming(Window* w)
{
    const unsigned char num_of_coef = 2;
    const double alpha = 0.54, /* also encountered as 0.54, 25/46, 0.54576, 0.53836, 0.53856 */
    			coef[2] = {alpha, 1-alpha}; /* see [1-3] for more info */

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL14                                                                                */
/********************************************************************************************/
void nuttall14(Window* w)
{
    const unsigned char num_of_coef = 3;
    const double coef[3] = {0.4243801, 0.4973406, 0.0782793};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* NUTTALL15                                                                                */
/********************************************************************************************/
void nuttall15(Window* w)
{
    const unsigned char num_of_coef = 4;
    const double coef[4] = {0.3635819, 0.4891775, 0.1365995, 0.0106411};

    raised_cosine(w, coef, num_of_coef);
}


/********************************************************************************************/
/* RAISED_COSINE                                                                            */
/********************************************************************************************/
/* TODO: consider using simmetry to optimize this (don't forget that the derivative change signs) */
void raised_cosine(Window* w, const double coef[], const unsigned char num_of_coef)
{
    const double length = w->length;
    double t, win_aux, deriv_aux;
    unsigned int i;
	unsigned char k;

    for (i = 0, t = length/-2, w->total_area = 0;
    	 i < length;
    	 i++, t += 1)
    {
        for (k = 0, win_aux = deriv_aux = 0; k < num_of_coef; k++)
        {
        	/* Firstly calculates the window */

            win_aux += (coef[k] * cos(TWO_PI * k * t / length));


            /* Then its derivative */

            deriv_aux += (k * coef[k] * sin(TWO_PI * k * t / length));
        }

        w->total_area += (w->content[i] = win_aux / length);

        w->derivative[i] = ( deriv_aux * (-TWO_PI / (length * length) ) );
    }


    /* Finally, the coeficients are copied */

    w->apodization_function.num_of_parameters = num_of_coef;
    w->apodization_function.parameters = (double*)emalloc(num_of_coef * sizeof(double));
    memmove(w->apodization_function.parameters, coef, num_of_coef * sizeof(double));
}

/********************************************************************************************/
/* OLD_RAISED_COSINE                                                                        */
/********************************************************************************************/
/* supposedly "symmetric", but no in even-simmetry, after Harris, 78 */
void old_raised_cosine(Window* w, const double coef[], const unsigned char num_of_coef)
{
    const double length = w->length;
    double t, win_aux, deriv_aux;
    unsigned int i;
	unsigned char k;

    for (i = 0, t = (length/-2 + 0.5), w->total_area = 0;
    	 i < length;
    	 i++, t += 1)
    {
        for (k = 0, win_aux = deriv_aux = 0; k < num_of_coef; k++)
        {
        	/* Firstly calculates the window */

            win_aux += (coef[k] * cos(TWO_PI * k * t / length));


            /* Then its derivative */

            deriv_aux += (k * coef[k] * sin(TWO_PI * k * t / length));
        }

        w->total_area += (w->content[i] = win_aux / length);

        w->derivative[i] = ( deriv_aux * (-TWO_PI / (length * length) ) );
    }


    /* Finally, the coeficients are copied */

    w->apodization_function.num_of_parameters = num_of_coef;
    w->apodization_function.parameters = (double*)emalloc(num_of_coef * sizeof(double));
    memmove(w->apodization_function.parameters, coef, num_of_coef * sizeof(double));
}

/********************************************************************************************/
/* GAUSSIAN                                                                                 */
/********************************************************************************************/
void gaussian(Window* w)
{
	gaussian_param(w, DEFAULT_GAUSSIAN_WINDOW_ALPHA);
}

/********************************************************************************************/
/* GAUSSIAN_PARAM                                                                        */
/********************************************************************************************/
void gaussian_param(Window* w, double alpha)
{
    const double length = w->length, ratio = -4*alpha/(length*length);
	unsigned int i;
	double t;

    for (i = 0, t = (length/-2 + 0.5), w->total_area = 0;
    	 i < length;
    	 i++, t += 1)
    {
	 	w->content[i] = exp(t*t * ratio);
	 	printf("win: %g\n", w->content[i]);
	 	w->total_area += w->content[i];

	 	w->derivative[i] = 2*ratio*exp(t*t * ratio);
	 	printf("deriv: %g\n", w->derivative[i]);
	}

    /* Finally, the coeficients are copied */

    w->apodization_function.num_of_parameters = 1;
    w->apodization_function.parameters = (double*)emalloc(sizeof(double));
    w->apodization_function.parameters[0] = alpha;
}

/********************************************************************************************/
/* HANNING_POISSON                                                                       */
/********************************************************************************************/
void hanning_poisson(Window* w)
{
	hanning_poisson_param(w, DEFAULT_HANNING_POISSON_WINDOW_ALPHA)	;
}

/********************************************************************************************/
/* HANNING_POISSON_PARAM                                                                 */
/********************************************************************************************/
void hanning_poisson_param(Window* w, double alpha)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;

	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
	 	w->content[i] = 0.5*(1+cos(PI*n/(win_length/2)))
							* exp(-alpha*fabs((double)n)/(win_length/2));
	 	w->total_area += w->content[i];
	}
}

/********************************************************************************************/
/* HELIE_A_W1                                                                            */
/********************************************************************************************/
void helie_a_w1(Window* w)
{
	helie_a_param(w, pow(10, -4), 21.6);
}

/********************************************************************************************/
/* HELIE_A_W6                                                                            */
/********************************************************************************************/
void helie_a_w6(Window* w)
{
	helie_a_param(w, 1.8, 0.92);
}

/********************************************************************************************/
/* HELIE_A_PARAM                                                                         */
/********************************************************************************************/
void helie_a_param(Window* w, double a, double b)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;
	double aux;


	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
 		aux = fabs(n)/(win_length/2);
	 	w->content[i] = pow(1 - aux, a)*exp(-b*pow(aux, 2));
	 	w->total_area += w->content[i];
	}
}


/********************************************************************************************/
/* Modified Barlett-Hann                                                                    */
/********************************************************************************************/
void mod_barlett_hann(Window* w)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;
	double t;


	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
 		t = n/(float)win_length;
	 	w->content[i] = 0.62 - 0.48 * fabs(t) + 0.38*cos(TWO_PI*t);
	 	w->total_area += w->content[i];
	}
}

/********************************************************************************************/
/* Bohman                                                                                   */
/********************************************************************************************/
void bohman(Window* w)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;
	double aux;


	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
 		aux = 2*fabs(n)/win_length;
	 	w->content[i] = (1-aux) * cos(PI*aux) + 1/PI * sin(PI*aux);
	 	w->total_area += w->content[i];
	}
}


/********************************************************************************************/
/* Bisquare                                                                                 */
/********************************************************************************************/
/* TODO: make it callable */
/* [1] http://www.clecom.co.uk/science/autosignal/help/Data_Tapering_Windows.htm */
void bisquare(Window* w)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;
	double t;


	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
 		t = n/(float)win_length;
	 	w->content[i] = pow(1.0-pow(fabs(i-0.5*win_length+0.5),2)/pow((0.5*win_length-0.5),2),2);
	 	w->total_area += w->content[i];
	}
}

/********************************************************************************************/
/* Welch                                                                                    */
/********************************************************************************************/
/* TODO: make it callable */
/* [1] http://www.clecom.co.uk/science/autosignal/help/Data_Tapering_Windows.htm */
void welch(Window* w)
{
	const unsigned int win_length = w->length;

	signed int n;
	unsigned int i;
	double t;


	for(i = 0, n = win_length/2 * (-1), w->total_area = 0; i < win_length; i++, n++)
 	{
 		t = n/(float)win_length;
	 	w->content[i] = 1.0-(((win_length-1)-2*i)/(win_length-1)*((win_length-1)-2*i)/(win_length-1));
	 	w->total_area += w->content[i];
	}
}

/********************************************************************************************/
/* OLD_HANN_WIN_MAG_FREQ_RESP                                                               */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
double old_hann_win_mag_freq_resp(double f) /* INCORRECT */
{
    if (f == 0) return 0.5;

	const double w = f * TWO_PI, sin_half_w = sin(w*0.5);

	/* the formula was analytically derived by the author
	 * and then computationally optimized by him */
	return 2*sin_half_w/w - sin_half_w/(TWO_PI - w) + sin_half_w/(TWO_PI + w);
}


/********************************************************************************************/
/* HANN_WIN_MAG_FREQ_RESP                                                                   */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
double hann_win_mag_freq_resp(double f)
{
    if (f == 0) return 0.5;

    const double fpow2 = f*f;

	/* the formula was analytically derived by the author
	 * and then computationally optimized by him */
	return f/PI*sin(PI*f)*(0.5/fpow2 - 0.5/(fpow2-1));
}


/********************************************************************************************/
/* HANN_WIN_LOG_POWER_FREQ_RESP                                                             */
/********************************************************************************************/
double hann_win_log_power_freq_resp(double f)
{
    double aux = hann_win_mag_freq_resp(f);

	return 10 * log10(aux*aux);
}


/********************************************************************************************/
/* RECTANGULAR_WIN_MAG_FREQ_RESP                                                            */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
double rectangular_win_mag_freq_resp(double f)
{
    if (f == 0) return 0.5;

	const double w = f * TWO_PI;

	return sin(w)/w;
}



/********************************************************************************************/
/* RAISED_COSINE_MAG_FREQ_RESP (based on Nuttall's paper)                                   */
/********************************************************************************************/
/* expect freq in bins */
double raised_cosine_mag_resp(const double lf, Window *window)
{
    if (lf == 0) return window->apodization_function.parameters[0];

    unsigned char k;
	double sign, sum;
	const double l2f2 = lf * lf;

    for (k = 0, sign = 1, sum = 0;
    	 k < window->apodization_function.num_of_parameters;
    	 k++, sign *= -1)
    {
	   	sum += (sign * window->apodization_function.parameters[k] / (l2f2 - k*k) );
    }

	return (lf / PI * sin(PI * lf) * sum);
}


/********************************************************************************************/
/* HANN_DERIVATIVE                                                                          */
/********************************************************************************************/
void hann_derivative(Window* w)
{
    const unsigned int length = w->length, last_index = length - 1;

    unsigned int i;
    double t;

    for (i = 0, t = (-1) * ((double)last_index/2), w->total_area = 0;
    	 i < length;
    	 i++, t += 1)
    {
        w->total_area += (w->content[i] = -sin(TWO_PI*t/last_index)*PI/last_index);
    }
}


/********************************************************************************************/
/* WINDOW_MAG_RESP                                                                          */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
double window_mag_resp(const double freq, Window* window)
{
	if (window->apodization_function.type == RECTANGULAR) { /* this shouldn't be necessary */
		if (freq == 0) return 1;      /* but raised cosine seems not to work for this case */

		double norm_freq = freq * INVERSE_OF_PI;

		return sin(norm_freq)/norm_freq;
	}
	if (window->apodization_function.type > RECTANGULAR
		&& window->apodization_function.type <= NUTTALL15)
	{
		return ( raised_cosine_mag_resp(freq, window) / raised_cosine_mag_resp(0, window) );
	}

	return 1;
}


/********************************************************************************************/
/* WINDOW_LOG_POWER_RESP                                                                    */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
double window_log_power_resp(const double freq, Window* window)
{
    double aux = window_mag_resp(freq, window)/window_mag_resp(0, window);

	return 10 * log10(aux*aux);
}


/********************************************************************************************/
/* INIT_UNPACK_SIN_TABLE                                                                    */
/********************************************************************************************/
void init_unpack_sin_table(const unsigned int n) {
	const unsigned int half_n = n/2;
    unsigned int k;

    unpack_sin_table = (double*)emalloc((half_n) * sizeof(double));

    for (k = 0; k < half_n; k++)
        unpack_sin_table[k] = sin(TWO_PI * k / n);
}


/********************************************************************************************/
/* FREE_UNPACK_SIN_TABLE                                                                    */
/********************************************************************************************/
void free_unpack_sin_table(void) {
	if (unpack_sin_table) free(unpack_sin_table);
}


/********************************************************************************************/
/* INIT_UNPACK_COS_TABLE                                                                    */
/********************************************************************************************/
void init_unpack_cos_table(const unsigned int n) {
	const unsigned int half_n = n/2;
    unsigned int k;

    unpack_cos_table = (double*)emalloc((half_n) * sizeof(double));

    for (k = 0; k < half_n; k++)
        unpack_cos_table[k] = cos(TWO_PI * k / n);
}


/********************************************************************************************/
/* FREE_UNPACK_COS_TABLE                                                                    */
/********************************************************************************************/
void free_unpack_cos_table(void) {
	if (unpack_cos_table) free(unpack_cos_table);
}


/********************************************************************************************/
/* PACK_REAL_IN_COMPLEX                                                                     */
/********************************************************************************************/
/* packs a real array into a half-sized complex array before the DFT */
void pack_real_in_complex(Complex packed[], const double unpacked[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int i, j, k;

	for (i = 0, j = 0, k = 1; i < half_n; i++, j+=2, k+=2) {
		packed[i].re = unpacked[j];
		packed[i].im = unpacked[k];
	}
}


/********************************************************************************************/
/* CLASSIC_PACK_REAL_IN_COMPLEX                                                             */
/********************************************************************************************/
/* packs a real array into a half-sized complex array before the DFT */
void classic_pack_real_in_complex(Complex packed[],
									const double unpacked[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int i;

	for (i = 0; i < half_n; i++) {
		packed[i].re = unpacked[2*i];
		packed[i].im = unpacked[2*i + 1];
	}
}

/********************************************************************************************/
/* REFERENCE_UNPACK_REAL_FROM_COMPLEX                                                       */
/********************************************************************************************/
/* unpacks a complex array into a half-sized one after the DFT */
void reference_unpack_real_from_complex(Complex unpacked[],
								const Complex packed[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int k;

	for (k = 1; k < half_n; k++) {
		c_copy(
			c_subtract(
				c_mult_by_real(c_add(packed[k], c_conjugate(packed[half_n - k])), 0.5),
				c_mult(
					c_mult(c_subtract(packed[k], c_conjugate(packed[half_n - k])),
						c_imaginary_to_complex(0.5)),
					c_exp(TWO_PI*k/n)
				)
			),
			&unpacked[k-1]);
	}
	c_copy(
		c_subtract(
			c_mult_by_real(c_add(packed[0], c_conjugate(packed[0])), 0.5),
			c_mult(
				c_mult(c_subtract(packed[0], c_conjugate(packed[0])),
					c_imaginary_to_complex(0.5)),
				c_exp(TWO_PI*half_n/n)
			)
		),
		&unpacked[half_n-1]);
}


/********************************************************************************************/
/* UNPACK_REAL_FROM_COMPLEX                                                                 */
/********************************************************************************************/
/* unpacks a complex array into a half-sized one after the DFT */
void unpack_real_from_complex(Complex unpacked[],
								const Complex packed[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int k;

	for (k = 1; k < half_n; k++) {
		const double cosine = unpack_cos_table[k], sine = unpack_sin_table[k];
		const unsigned int half_n_minus_k = half_n - k;

		unpacked[k-1].re = 0.5 * (
								packed[k].re * (1 + sine)
								+
								packed[half_n_minus_k].re * (1 - sine)
								+
								cosine * (packed[k].im + packed[half_n_minus_k].im)
								);
		unpacked[k-1].im = 0.5 * (
								packed[k].im * (1 + sine)
								+
								packed[half_n_minus_k].im * (sine - 1)
								-
								cosine * (packed[k].re - packed[half_n_minus_k].re)
								);
	}
	unpacked[half_n-1].re = packed[0].re - packed[0].im;
	unpacked[half_n-1].im = 0;
}


/********************************************************************************************/
/* CLASSIC_UNPACK_REAL_FROM_COMPLEX                                                         */
/********************************************************************************************/
/* unpacks a complex array into a half-sized one after the DFT */
void classic_unpack_real_from_complex(Complex unpacked[],
								const Complex packed[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int k;

	for (k = 1; k < half_n; k++) {
		unpacked[k-1].re = 0.5 * (
								packed[k].re * (1 + unpack_sin_table[k])
								+
								packed[half_n - k].re * (1 - unpack_sin_table[k])
								+
								unpack_cos_table[k] * (packed[k].im + packed[half_n - k].im)
								);
		unpacked[k-1].im = 0.5 * (
								packed[k].im * (1 + unpack_sin_table[k])
								+
								packed[half_n - k].im * (unpack_sin_table[k] -1)
								-
								unpack_cos_table[k] * (packed[k].re - packed[half_n - k].re)
								);
	}

	unpacked[half_n-1].re = packed[0].re - packed[0].im;
	unpacked[half_n-1].im = 0;
}


/********************************************************************************************/
/* PREPARE_FFT                                                                              */
/********************************************************************************************/
void prepare_fft(const unsigned int n)
{
	init_rev_table(n);
	init_sin_table(n);
	init_cos_table(n);
}

/********************************************************************************************/
/* FINALIZE_FFT                                                                             */
/********************************************************************************************/
void finalize_fft(void)
{
    free_rev_table();
    free_sin_table();
    free_cos_table();
}

/********************************************************************************************/
/* PREPARE_FFT_REAL_SIGNAL_TO_POWER_PHASE                                                   */
/********************************************************************************************/
void prepare_fft_real_signal_to_power_phase(const unsigned int n)
{
	init_unpack_cos_table(n);
	init_unpack_sin_table(n);

	packed = (Complex*)emalloc(n/2 * sizeof(Complex));
	unpacked = (Complex*)emalloc(n/2 * sizeof(Complex));

	prepare_fft(n/2);
}


/********************************************************************************************/
/* FINALIZE_FFT_REAL_SIGNAL_TO_POWER_PHASE                                                  */
/********************************************************************************************/
void finalize_fft_real_signal_to_power_phase(void)
{
	free_unpack_cos_table();
	free_unpack_sin_table();

	if (packed) free(packed);
	if (unpacked) free(unpacked);

	finalize_fft();
}

/********************************************************************************************/
/* PREPARE_FFT_REAL_SIGNAL                                                                  */
/********************************************************************************************/
void prepare_fft_real_signal(const unsigned int n)
{
	init_unpack_cos_table(n);
	init_unpack_sin_table(n);

	packed = (Complex*)emalloc(n/2 * sizeof(Complex));
	unpacked = (Complex*)emalloc(n/2 * sizeof(Complex));

	prepare_fft(n/2);
}


/********************************************************************************************/
/* FINALIZE_FFT_REAL_SIGNAL                                                                 */
/********************************************************************************************/
void finalize_fft_real_signal(void)
{
	free_unpack_cos_table();
	free_unpack_sin_table();

	if (packed) free(packed);
	if (unpacked) free(unpacked);

	finalize_fft();
}


/********************************************************************************************/
/* FFT_REAL_SIGNAL_TO_POWER_PHASE                                                           */
/********************************************************************************************/
void fft_real_signal_to_power_phase(double power[], double phase[],
									const double signal[], const unsigned int n)
{
	const unsigned int half_n = n/2;
	unsigned int i;

	pack_real_in_complex(unpacked, signal, n);

	fft(unpacked, packed, half_n);

	unpack_real_from_complex(unpacked, packed, n);

	if (power != NULL) {
		if (phase != NULL) {
			for (i = 0; i < half_n; i++) {
				power[i] = unpacked[i].re*unpacked[i].re + unpacked[i].im*unpacked[i].im;
				phase[i] = atan2(unpacked[i].im, unpacked[i].re);
			}
			power[half_n-1] /= 4;
		} else {
			for (i = 0; i < half_n; i++) {
				power[i] = unpacked[i].re*unpacked[i].re + unpacked[i].im*unpacked[i].im;
			}
			power[half_n-1] /= 4;
		}
	} else {
		if (phase != NULL) {
			for (i = 0; i < half_n; i++) {
				phase[i] = atan2(unpacked[i].im, unpacked[i].re);
			}
		}
	}
}


/********************************************************************************************/
/* FFT_REAL_SIGNAL                                                                          */
/********************************************************************************************/
void fft_real_signal(Complex spectrum[], const double signal[], const unsigned int n)
{
	const unsigned int half_n = n/2;

	pack_real_in_complex(spectrum, signal, n);

	fft(spectrum, packed, half_n);

	unpack_real_from_complex(spectrum, packed, n);
}


/********************************************************************************************/
/* REFERENCE_FFT_REAL_SIGNAL_POWER_PHASE                                                    */
/********************************************************************************************/
void reference_fft_real_signal_power_phase(const float s[], const unsigned int m,
				const Apodization_Function_Type apod, const unsigned int zeroes,
				double pow[], double ph[])
{
	const unsigned int n = m + zeroes;
	if (!is_power_of_2(n)) {
		fprintf(stderr, "Error in reference_fft_real_signal_power_phase():\n"
							"\tframe total size must be a power of two.\n");
		return;
	}

	if (pow == NULL && ph == NULL) return;
	/*
	const unsigned int half_n = n/2;
	Complex *aux = (Complex*)emalloc(n/2*sizeof(Complex)),
			*packed_spectrum = (Complex*)emalloc(n/2*sizeof(Complex));

	pack_real_in_complex(aux, s, n);

	reference_fft(aux, packed_spectrum, half_n);

	reference_unpack_real_from_complex(aux, packed_spectrum, n);
	*/
	Complex *sig, *spec;
	Window w;

	w.length = m;
	w.apodization_function.type = apod;
	init_window(&w);

	sig = (Complex*)emalloc(n*sizeof(Complex));
	spec = (Complex*)emalloc(n*sizeof(Complex));

	unsigned int i;
	for (i = 0; i < m; i++) {
		/*
		sig[i].re = s[i] * w.content[i];
		*/
		sig[i].re = s[i];
		sig[i].im = 0;
	}
	for (i = m; i < n; i++) {
		sig[i].re = sig[i].im = 0;
	}
	reference_fft(sig, spec, n);

	if (pow) {
		for (i = 0; i < n/2; i++) {
			pow[i] = c_modulus_square(spec[i+n/2-1]);
		}
	}
	if (ph) {
		for (i = 0; i < n/2; i++) {
			ph[i] = c_phase(spec[i+n/2-1]);
		}
	}
	free(sig);
	free(spec);
	free_window(&w);

	/*
	if (pow) {
		if (ph) {
			for (i = 0; i < half_n; i++) {
				pow[i] = c_modulus_square(aux[i]);
				ph[i] = c_phase(aux[i]);
			}
		} else {
			for (i = 0; i < half_n; i++) {
				pow[i] = c_modulus_square(aux[i]);
			}
		}
	} else {
		for (i = 0; i < half_n; i++) {
			ph[i] = c_phase(aux[i]);
		}
	}
	*/
}
