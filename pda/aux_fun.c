#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>

#include "numeric_constants.h"
#include "data_structures.h"


/********************************************************************************************/
/* BOUND                                                                                    */
/********************************************************************************************/
void bound(unsigned char *var, unsigned char min_val, unsigned char max_val)
{
	if (*var < min_val) *var = min_val;
	if (*var > max_val) *var = max_val;
}

/********************************************************************************************/
/* FBOUND                                                                                   */
/********************************************************************************************/
void fbound(double *var, double min_val, double max_val)
{
	if (*var < min_val) *var = min_val;
	if (*var > max_val) *var = max_val;
}

/********************************************************************************************/
/* PARSE_TRUE_OR_FALSE                                                                      */
/********************************************************************************************/
signed char parse_true_or_false(const char argument[])
{
		if (!strcmp(argument, "true")) {
			return TRUE;
		} else if (!strcmp(argument, "false")) {
			return FALSE;
		} else {
			return ERROR;
		}
}


/********************************************************************************************/
/* VERIFY_NONNEGATIVENESS                                                                   */
/********************************************************************************************/
void verify_nonnegativeness(const double value, const char varname[], signed char *status)
{
	if (value < 0) {
        	fprintf(stderr, "Error: %s must be non-negative.\n", varname);
        	*status = ERROR;
	}
}


/********************************************************************************************/
/* VERIFY_POSITIVENESS                                                                      */
/********************************************************************************************/
void verify_positiveness(const double value, const char varname[], signed char *status)
{
	if (value <= 0) {
        	fprintf(stderr, "Error: %s must be strictly positive.\n", varname);
        	*status = ERROR;
	}
}


/********************************************************************************************/
/* VERIFY_NOT_SMALLER_THAN                                                                  */
/********************************************************************************************/
void verify_not_smaller_than(const double var1, const double var2, const char varname1[],
												const char varname2[], signed char *status)
{
	if (var1 < var2) {
        	fprintf(stderr, "Error: %s must be equal or greater than %s.\n", varname1, varname2);
        	*status = ERROR;
	}
}


/********************************************************************************************/
/* PARSE_ARGUMENTS                                                                          */
/********************************************************************************************/
signed char parse_arguments(const int n_args, char** argument, FILE **F0_fp, FILE **onset_fp)
{
    signed char status = OK;
	int i;

	/* set default values */	

	/* parse arguments */

	for (i = 1; i < n_args; i++)
	{
		if (!strcmp(argument[i], "--F0_input_file") || !strcmp(argument[i], "-f"))
		{
      		if (++i < n_args)
        	{
			    if ((*F0_fp = fopen(argument[i], "r")) == NULL)
			    {
		        fprintf(stderr, "Error opening F0 input file '%s'.\n", argument[i]);
		        status = ERROR;
			    }
        	}
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else
  		{
        	fprintf(stderr, "Unrecognized parameter: '%s'\n", argument[i]);
        	status = ERROR;
        }
	}
	
	/* validate data */
	
	return status;
}

/********************************************************************************************/
/* EMALLOC                                                                                  */
/********************************************************************************************/
void* emalloc(const size_t size)
{
    void* p;
    
    if ( ( p = malloc(size) ) != NULL )
    	return p;
    
    fprintf(stderr, "Error allocating memory. Aborting...\n");
    exit ( ERROR );
}

/********************************************************************************************/
/* ECALLOC                                                                                  */
/********************************************************************************************/
void* ecalloc(const size_t size)
{
    void* p;

    if ( ( p = calloc(1, size) ) != NULL )
    	return p;

    fprintf(stderr, "Error allocating memory. Aborting...\n");
    exit ( ERROR );
}


/********************************************************************************************/
/* EREALLOC                                                                                 */
/********************************************************************************************/
void* erealloc(void *p, const size_t size)
{
    void* q = realloc(p, size);
    
    if ( q != NULL )
    	return q;
    	
    fprintf(stderr, "Error reallocating memory. Aborting...\n");
    exit ( ERROR );
}


/********************************************************************************************/
/* LOG2                                                                                     */
/********************************************************************************************/
double log2(double x)
{
	return ( log10(x) / LOGARITHM_OF_2_IN_BASE_10 );
}


/********************************************************************************************/
/* FREE_2D                                                                                  */
/********************************************************************************************/
void free_2d(void **p, const unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		free(p[i]);
	free(p);
}

/********************************************************************************************/
/* HERTZ_TO_BARK                                                                            */
/********************************************************************************************/
/* performs Hertz frequency conversion to Bark Scale
 * according to Traunmüller (1983, 1988, 1990) */
double hertz_to_bark(const double f)
{
	double bark = (26.81*f)/(1960+f)-0.53;
	
	if (bark < 2) {
		bark = bark + 0.15*(2-bark);
	} else if (bark > 20.1) {
		bark = bark + 0.22*(bark-20.1);
	}
	
	return bark;
}

/********************************************************************************************/
/* HERTZ_TO_ERB                                                                             */
/********************************************************************************************/
/* performs Hertz frequency conversion to Bark Scale
 * according to Traunmüller (1983, 1988, 1990) */
double hertz_to_erb(const double f)
{
	return 21.4 * log10(4.37e-3*f + 1);
}
