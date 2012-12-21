#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "filters.h"
#include "config_variables.h"
#include "config_constants.h"
#include "numeric_constants.h"
#include "aux_fun.h"


int main(int argc, char *argv[])
{
	unsigned int num_of_peaks = 0, i;
  	double last_time = -1, curr_time, curr_freq, curr_amp, *freq, *amp;
  	char *line = (char*) emalloc(MAX_LINE_LENGTH * sizeof(char));

  	FILE *fp;
  
  	if (parse_arguments(argc, argv, &fp) != OK) {
  		return ERROR;
  	}

  	freq = (double*)emalloc(MAX_NUM_OF_PEAKS * sizeof(double));
  	amp = (double*)emalloc(MAX_NUM_OF_PEAKS * sizeof(double));
  
	while (fgets(line, MAX_LINE_LENGTH, fp)) {
		if (line[0] == COMMENT_CHARACTER) continue;
  
		if (sscanf(line, "%lf%lf%lf", &curr_time, &curr_freq, &curr_amp) != 3) continue;  
  
		if (curr_time != last_time) {
			
		    if (num_of_peaks != 0)
	   		{
	    		num_of_peaks = multi_filter(amp, freq, num_of_peaks,
	    								MIN_FREQUENCY, MAX_FREQUENCY, MIN_ABSOLUTE_MAGNITUDE,
	    								MAX_MAGNITUDE_GAP, MAX_ADJACENT_MAGNITUDE_GAP);
	
		    	for (i = 0; i < num_of_peaks; i++) {
		    		printf("%g %g %g\n", last_time, freq[i], amp[i]);
	        	}
	    	}
	  		last_time = curr_time;
	  		num_of_peaks = 0;
		}
		freq[num_of_peaks] = curr_freq;
		amp[num_of_peaks] = curr_amp;
		
		num_of_peaks++;
	}
	fclose(fp);
   
   
	if (num_of_peaks != 0)
	{
		num_of_peaks = multi_filter(amp, freq, num_of_peaks,
								MIN_FREQUENCY, MAX_FREQUENCY, MIN_ABSOLUTE_MAGNITUDE,
								MAX_MAGNITUDE_GAP, MAX_ADJACENT_MAGNITUDE_GAP);

		for (i = 0; i < num_of_peaks; i++) {
			printf("%g %g %g\n", last_time, freq[i], amp[i]);
		}
  	}

	free(freq);
	free(amp);
  
  	return OK;
}
