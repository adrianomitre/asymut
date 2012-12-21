#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "numeric_constants.h"
#include "config_variables.h"
#include "config_constants.h"
#include "numeric_constants.h"


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
signed char parse_arguments(const int n_args, char **argument, FILE **fp)
{
	unsigned int i;
    signed char status = OK;
	
	/* set default values */	
	
	*fp = stdin;
	
	MIN_FREQUENCY = DEFAULT_MIN_FREQUENCY;
	MAX_FREQUENCY = DEFAULT_MAX_FREQUENCY;
	
	MIN_ABSOLUTE_MAGNITUDE = DEFAULT_MIN_ABSOLUTE_MAGNITUDE;
	MAX_MAGNITUDE_GAP = DEFAULT_MAX_MAGNITUDE_GAP;
	MAX_ADJACENT_MAGNITUDE_GAP = DEFAULT_MAX_ADJACENT_MAGNITUDE_GAP;
	
	/* parse arguments */

	for (i = 1; i < n_args; i++)
	{
		if (!strcmp(argument[i], "--input_file") || !strcmp(argument[i], "-i"))
		{
      		if (++i < n_args)
        	{
			    if ((*fp = fopen(argument[i], "r")) == NULL)
			    {
		        fprintf(stderr, "Error opening input file '%s'.\n", argument[i]);
		        status = ERROR;
			    }
        	}
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_freq") || !strcmp(argument[i], "-l"))
		{
      		if (++i < n_args) MIN_FREQUENCY = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_freq") || !strcmp(argument[i], "-h"))
		{
      		if (++i < n_args) MAX_FREQUENCY = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_mag") || !strcmp(argument[i], "-m"))
		{
      		if (++i < n_args) MIN_ABSOLUTE_MAGNITUDE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_gap") || !strcmp(argument[i], "-g"))
		{
      		if (++i < n_args) MAX_MAGNITUDE_GAP = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_adj_gap") || !strcmp(argument[i], "-a"))
		{
      		if (++i < n_args) MAX_ADJACENT_MAGNITUDE_GAP = atof(argument[i]);
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

	verify_positiveness(MIN_FREQUENCY, "minimum frequency", &status);	
	verify_positiveness(MAX_FREQUENCY, "maximum frequency", &status);
	verify_not_smaller_than(MAX_FREQUENCY, MIN_FREQUENCY,
										"maximum frequency", "minimum frequency", &status);

	verify_nonnegativeness(MAX_MAGNITUDE_GAP, "maximum magnitude gap", &status);
	verify_nonnegativeness(MAX_ADJACENT_MAGNITUDE_GAP, "maximum adjacent magnitude gap", &status);
	
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
