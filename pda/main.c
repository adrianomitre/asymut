#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "f0_estimation_methods.h"
#include "numeric_constants.h"
#include "config_constants.h"
#include "aux_fun.h"

/***********************************************************************
 * MAIN FUNCTION                                                       *
 ***********************************************************************/
int main(int argc, char *argv[])
{
  double last_time = -1, curr_time, curr_freq, curr_amp, *freq, *amp,
     	min_f0 = DEFAULT_MIN_F0,
      	max_f0 = DEFAULT_MAX_F0,
      	minRelativeEnergy = DEFAULT_MIN_RELATIVE_ENERGY,
      	maxRelativeEnergyAbove = DEFAULT_MAX_RELATIVE_ENERGY_ABOVE,
      	maxRelativeEnergyBelow = DEFAULT_MAX_RELATIVE_ENERGY_BELOW,
       	f0_est;
  unsigned int num_of_peaks = 0;
  char *line = (char*) emalloc(MAX_LINE_LENGTH * sizeof(char));
  FILE *fp;
  EstimationParameters *ep = (EstimationParameters*) emalloc(sizeof(EstimationParameters));
  Boolean eof_reached;

  if (argc < 2) {
    /*
      fprintf(stderr, "Error: user must provide input file name.\n");
      fprintf(stderr, "Usage: '%s' <filename> [min_f0] [max_f0] [MRP] [MRP_below] [MRP_above].\n", argv[0]);
      return 1;
    */
  	fp = stdin;
  } else {
    fp = fopen(argv[1], "r");
    if (!fp) {
      fprintf(stderr, "Error opening file '%s'\n", argv[1]);
      return 1;
    }
  }

  if (argc >= 3) {
      min_f0 = atof(argv[2]);
      if (min_f0 <= 0) {
          fprintf(stderr, "Error: min_f0 must be positive.\n");
          return 1;
      }
  }
  if (argc >= 4) {
      max_f0 = atof(argv[3]);
      if (max_f0 <= 0) {
          fprintf(stderr, "Error: max_f0 must be positive.\n");
          return 1;
      }
  }
  if (argc >= 5) {
      minRelativeEnergy = atof(argv[4]);
      if (minRelativeEnergy < 0 || minRelativeEnergy > 1) {
          fprintf(stderr, "Error: minRelativeEnergy must be within 0 and 1.\n");
          return 1;
      }
  }
  if (argc >= 6) {
      maxRelativeEnergyBelow = atof(argv[5]);
      if (maxRelativeEnergyBelow < 0 || maxRelativeEnergyBelow > 1) {
          fprintf(stderr, "Error: maxRelativeEnergyBelow must be within 0 and 1.\n");
          return 1;
      }
  }
  if (argc >= 7) {
      maxRelativeEnergyAbove = atof(argv[6]);
      if (maxRelativeEnergyAbove < 0 || maxRelativeEnergyAbove > 1) {
          fprintf(stderr, "Error: maxRelativeEnergyAbove must be within 0 and 1.\n");
          return 1;
      }
  }


  freq = (double*) emalloc(MAX_PEAKS * sizeof(double));
  amp = (double*) emalloc(MAX_PEAKS * sizeof(double));

  /*
  if (!strcmp(argv[1], "stdin")) {
  	fp = stdin;
  }
  else {
      if (!(fp = fopen(argv[1], "r"))) {
      	fprintf(stderr, "Error opening file '%s'\n", argv[1]);
        return 1;
      }
  }
  */

  prepare_critical_band_sum(MAX_PEAKS);

  ep->minAbsFreq = min_f0;
  ep->maxAbsFreq = max_f0;
  ep->maxRelativeEnergyBelow = maxRelativeEnergyBelow;
  ep->maxRelativeEnergyAbove = maxRelativeEnergyAbove;
  ep->minRelativeEnergy = minRelativeEnergy;
  ep->maxRelativeInharmonicity = 0.5*(HALFTONE_UP - 1);
  ep->opMode = PRINT_REFINED_F0_ESTIMATE;

  eof_reached = FALSE;
  while(!eof_reached)
  {
    if (fgets(line, MAX_LINE_LENGTH, fp) == NULL) {
      curr_time = DBL_MAX;
      eof_reached = TRUE;
    } else {
      if (line[0] == COMMENT_CHARACTER) continue;

      if (sscanf(line, "%lf%lf%lf", &curr_time, &curr_freq,
              &curr_amp) != 3) continue;
    }

		if (curr_time != last_time) {

	        if (num_of_peaks != 0) {
    	        f0_est = mitre_f0_est(freq, amp, num_of_peaks, last_time, ep);

	  			/*
				f0_est = full_average_f0_est(freq, amp, num_of_peaks, min_f0, max_f0, 0, 0);

	            if (f0_est > 0) {
	            	printf("%g %g\n", last_time, f0_est);
	            } else {
	            	fprintf(stderr, "%g %g\n", last_time, f0_est);
	            }
	            */

	        }

      	last_time = curr_time;
      	num_of_peaks = 0;
      }
	if (curr_amp > 0) { /* No need for filtering */
		freq[num_of_peaks] = curr_freq;
		amp[num_of_peaks] = curr_amp;

		num_of_peaks++;
	}
  }
  fclose(fp);

  finalize_critical_band_sum();

  return 0;
}

