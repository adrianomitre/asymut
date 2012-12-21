#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>

#include "numeric_constants.h"
#include "config_constants.h"
#include "config_variables.h"


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

	MIN_NOTE_DURATION = DEFAULT_MIN_NOTE_DURATION;
	MAX_DELAY_AFTER_ONSET = DEFAULT_MAX_DELAY_AFTER_ONSET;
	MAX_GAP_INSIDE_NOTE = DEFAULT_MAX_GAP_INSIDE_NOTE;
	
	MAX_ATTACK_LENGTH = DEFAULT_MAX_ATTACK_LENGTH;
	MIN_ATTACK_ABSOLUTE_INCREASE = DEFAULT_MIN_ATTACK_ABSOLUTE_INCREASE;
	MIN_ATTACK_RELATIVE_INCREASE = DEFAULT_MIN_ATTACK_RELATIVE_INCREASE;

	FREQ_REF_A4 = DEFAULT_FREQ_REF_A4;
	
	MAX_INTENSITY_LEVEL = DEFAULT_MAX_INTENSITY_LEVEL;
	MIN_RELATIVE_INTENSITY_LEVEL = DEFAULT_MIN_RELATIVE_INTENSITY_LEVEL;
	MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL = DEFAULT_MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL;
	
	PITCH_RANGE_LOWEST_NOTE = DEFAULT_PITCH_RANGE_LOWEST_NOTE;
	PITCH_RANGE_HIGHEST_NOTE = DEFAULT_PITCH_RANGE_HIGHEST_NOTE;

	PRINT_NOTE_AS = DEFAULT_NOTE_PRINTING_STYLE;
	
	TEMPO = DEFAULT_TEMPO_IN_BEATS_PER_MINUTE;
	INITIAL_PAUSE = DEFAULT_INITIAL_OFFSET_IN_SECONDS;
	INSTRUMENT = DEFAULT_GM_INSTRUMENT_NUMBER;
	CHANNEL_VOLUME = DEFAULT_GM_CHANNEL_VOLUME;
	
	ALLOW_NOTES_WITHOUT_ONSET = DEFAULT_NOTES_WITHOUT_ONSET_PERMISSION;
	
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
		else if (!strcmp(argument[i], "--onset_input_file") || !strcmp(argument[i], "-o"))
		{
      		if (++i < n_args)
        	{
			    if ((*onset_fp = fopen(argument[i], "r")) == NULL)
			    {
		        fprintf(stderr, "Error opening onset input file '%s'.\n", argument[i]);
		        status = ERROR;
			    }
        	}
      		else
        	{
                fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
                status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_duration") || !strcmp(argument[i], "-d"))
		{
      		if (++i < n_args) MIN_NOTE_DURATION = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_attack_length") || !strcmp(argument[i], "-k"))
		{
      		if (++i < n_args) MAX_ATTACK_LENGTH = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_attack_abs_incr"))
		{
      		if (++i < n_args) MIN_ATTACK_ABSOLUTE_INCREASE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_attack_rel_incr"))
		{
      		if (++i < n_args) MIN_ATTACK_RELATIVE_INCREASE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_delay") || !strcmp(argument[i], "-e"))
		{
      		if (++i < n_args) MAX_DELAY_AFTER_ONSET = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--max_gap") || !strcmp(argument[i], "-g"))
		{
      		if (++i < n_args) MAX_GAP_INSIDE_NOTE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--ref_a4_freq") || !strcmp(argument[i], "-r"))
        {
      		if (++i < n_args) FREQ_REF_A4 = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--lowest_note") || !strcmp(argument[i], "-l"))
        {
      		if (++i < n_args) PITCH_RANGE_LOWEST_NOTE = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--highest_note") || !strcmp(argument[i], "-h"))
        {
      		if (++i < n_args) PITCH_RANGE_HIGHEST_NOTE = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--print_mode") || !strcmp(argument[i], "-m"))
        {
      		if (++i < n_args)
      		{
            	if (!strcmp(argument[i], "list")) PRINT_NOTE_AS = HUMAN_READABLE_NOTE_LIST;
             	else if (!strcmp(argument[i], "midifiable-csv")) PRINT_NOTE_AS = MIDIFIABLE_CSV;
             	else if (!strcmp(argument[i], "piano-roll")) PRINT_NOTE_AS = GNUPLOT_ARROWS;
                else
                {
                    fprintf(stderr, "Unrecognized print mode: '%s'\n", argument[i]);
                    status = ERROR;
                }
      		}
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--tempo") || !strcmp(argument[i], "-t"))
        {
      		if (++i < n_args) TEMPO = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--initial_pause") || !strcmp(argument[i], "-p"))
        {
      		if (++i < n_args) INITIAL_PAUSE = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--instrument") || !strcmp(argument[i], "-i"))
        {
      		if (++i < n_args) INSTRUMENT = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--channel_volume") || !strcmp(argument[i], "-c"))
        {
      		if (++i < n_args) CHANNEL_VOLUME = atoi(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--allow_notes_without_onset") || !strcmp(argument[i], "-a"))
        {
      		if (++i < n_args) {
      			signed char aux = parse_true_or_false(argument[i]);
      			
      			if (aux == ERROR) {
             		fprintf(stderr, "Unrecognized true or false value: '%s'\n", argument[i]);
      				status = ERROR;
      			} else {
      				ALLOW_NOTES_WITHOUT_ONSET = (Boolean)aux;
      			}
      		}
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}		
        else if (!strcmp(argument[i], "--max_intensity") || !strcmp(argument[i], "-x"))
        {
      		if (++i < n_args) MAX_INTENSITY_LEVEL = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
        else if (!strcmp(argument[i], "--min_relative_intensity") || !strcmp(argument[i], "-v"))
        {
      		if (++i < n_args) MIN_RELATIVE_INTENSITY_LEVEL = atof(argument[i]);
      		else
        	{
             	fprintf(stderr, "Incomplete parameter: '%s'\n", argument[--i]);
             	status = ERROR;
            }
		}
		else if (!strcmp(argument[i], "--min_relative_restimulus_intensity") || !strcmp(argument[i], "-y"))
        {
      		if (++i < n_args) MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL = atof(argument[i]);
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
	
	verify_nonnegativeness(MIN_NOTE_DURATION, "minimum note duration", &status);
	verify_nonnegativeness(MAX_DELAY_AFTER_ONSET, "max delay after onset", &status);
	verify_nonnegativeness(MAX_GAP_INSIDE_NOTE, "max gap inside note", &status);

	verify_positiveness(FREQ_REF_A4, "A4 reference frequency", &status);

	verify_positiveness(TEMPO, "tempo", &status);
	verify_nonnegativeness(INITIAL_PAUSE, "initial pause", &status);
	
	verify_nonnegativeness(MAX_INTENSITY_LEVEL, "maximum intensity level", &status);
	
	bound(&INSTRUMENT, 0, 127);
	bound(&CHANNEL_VOLUME, 0, 127);
	bound(&PITCH_RANGE_LOWEST_NOTE, 0, 127);
	bound(&PITCH_RANGE_HIGHEST_NOTE, 0, 127);	
	fbound(&MIN_RELATIVE_INTENSITY_LEVEL, 0, 1);
	fbound(&MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL, 0, 1);

	verify_nonnegativeness(MAX_ATTACK_LENGTH, "maximum attack length", &status);
	verify_nonnegativeness(MIN_ATTACK_RELATIVE_INCREASE, "minimum attack relative increase", &status);
	verify_positiveness(MIN_ATTACK_ABSOLUTE_INCREASE, "minimum attack absolute increase", &status);
	
	verify_not_smaller_than(PITCH_RANGE_HIGHEST_NOTE, PITCH_RANGE_LOWEST_NOTE,
										"highest note", "lowest note", &status);
	
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
/* LOG2                                                                                     */
/********************************************************************************************/
double log2(double x)
{
	return ( log10(x) / LOGARITHM_OF_2_IN_BASE_10 );
}
