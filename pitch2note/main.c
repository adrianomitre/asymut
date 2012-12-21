#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "numeric_constants.h"
#include "config_constants.h"
#include "config_variables.h"
#include "aux_fun.h"
#include "notes.h"

/***********************************************************************
 * MAIN FUNCTION                                                       *
 ***********************************************************************/
int main(int argc, char **argv) {
	FILE *F0_fp = stdin, *onsets_fp = NULL;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * I. Parse arguments (input files)                                                      * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    if (parse_arguments(argc, argv, &F0_fp, &onsets_fp) != OK) {
    	return ERROR;
    }

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * II. Prepare for processing (initialize structures and parameters)                     *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	Pitch_Range *pr = NULL;
	
	if (init_pitch_range(&pr, PITCH_RANGE_LOWEST_NOTE, PITCH_RANGE_HIGHEST_NOTE) != OK)
		return ERROR;
		
	set_bpm(TEMPO);
	set_offset(INITIAL_PAUSE); /* initial 'pause' */
	set_instrument(INSTRUMENT);


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * III. Process files                                                                     * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	Flag flag = DOWN;
	double last_onset_time = DBL_MIN,
			next_onset_time = DBL_MAX,
			curr_F0_time = DBL_MIN,
			curr_F0_estimate = DBL_MIN,
			curr_F0_intensity = DBL_MIN;
	signed char curr_note_index ;
	char *line = (char*) emalloc(MAX_LINE_LENGTH * sizeof(char));
	

	if (onsets_fp != NULL)
	{
		/* discards comment lines at the beginning of onsets file */
		while (fgets(line, MAX_LINE_LENGTH, onsets_fp)) {
			if (line[0] == COMMENT_CHARACTER) continue;

			if (sscanf(line, "%lf", &next_onset_time) == 1) {
				break;
			}
		}
		if (next_onset_time == DBL_MAX && ALLOW_NOTES_WITHOUT_ONSET == FALSE) {
			fprintf(stderr, "Nothing to do: onsets file is empty.\n");
			return ERROR;
		}		
	}
	
	print_file_header();

	while (fgets(line, MAX_LINE_LENGTH, F0_fp))
	{
		if (line[0] == COMMENT_CHARACTER) continue;

		if (sscanf(line, "%lf%lf%lf", &curr_F0_time, &curr_F0_estimate,
						&curr_F0_intensity) != 3) continue;
		
		turn_off_obsolete_notes(pr, curr_F0_time);

		if ((curr_note_index = freq2note(curr_F0_estimate, pr)) < 0) continue;

		stimulate_note(&pr->note[curr_note_index], curr_F0_time, curr_F0_intensity);

		if (onsets_fp != NULL) {
			if (curr_F0_time >= next_onset_time) {
				last_onset_time = next_onset_time;
				
				if (!fgets(line, MAX_LINE_LENGTH, onsets_fp)
					|| sscanf(line, "%lf", &next_onset_time) != 1)
				{
					next_onset_time = DBL_MAX;
				}
				
				flag = UP;
			} else 	if ( (curr_F0_time - last_onset_time) > MAX_DELAY_AFTER_ONSET ) {
				flag = DOWN;
			}
		} else {
			last_onset_time = curr_F0_time;
		}
		
		if (flag == UP) {
			turn_on_new_notes(pr, last_onset_time, curr_F0_time);
		} else if (ALLOW_NOTES_WITHOUT_ONSET == TRUE) {
			turn_on_new_notes(pr, curr_F0_time, curr_F0_time);
		}
	}
	turn_off_obsolete_notes(pr, DBL_MAX);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * V. Print file footer                                                                 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	print_file_footer(curr_F0_time);
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * VI. Finalize (free allocated memory)                                                  *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free_pitch_range(&pr);


	return OK;
}
