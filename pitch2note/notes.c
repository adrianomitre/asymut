#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "notes.h"
#include "aux_fun.h"
#include "numeric_constants.h"
#include "config_constants.h"
#include "config_variables.h"

const static double TICS_PER_MIN = 60000000;
static unsigned int TICS_PER_BEAT;
static double TICS_PER_SEC;
static double TIME_OFFSET; /* in seconds */
static unsigned char GM_INSTRUMENT_NUMBER;
static double MIN_PRINTED_NOTE_EVIDENCE_LEVEL;
static double MAX_PRINTED_NOTE_EVIDENCE_LEVEL;
static signed char LOWEST_PRINTED_NOTE; /* lowest printed note NUMBER, actually */
static signed char HIGHEST_PRINTED_NOTE; /* highest printed note NUMBER, actually */

/********************************************************************************************/
/* INIT_PITCH_RANGE                                                                         */
/********************************************************************************************/
signed char init_pitch_range(Pitch_Range **pitch_range, const unsigned char first_note,
																const unsigned char last_note)
{
    unsigned char i, note_class;
    unsigned char const num_of_notes = last_note - first_note + 1;
    Pitch_Range *pr;

    if (last_note < first_note)
    {
        fprintf(stderr, "Error in Pitch_Range initialization: 'first_note' is greater than "
        																	"'last_note'.\n");
    	return ERROR;
   	}

	*pitch_range = (Pitch_Range*) emalloc(sizeof(Pitch_Range));

   	pr = *pitch_range;

    pr->note = (Note*) emalloc(num_of_notes * sizeof(Note));
    pr->num_of_notes = num_of_notes;
    pr->current_polyphony = 0;

    for (i = 0; i < num_of_notes; i++)
    {
        pr->note[i].fund_freq = FREQ_REF_A4 * pow(2, (double)(i + first_note - 69) / 12);
    	pr->note[i].last_evidence_time = DBL_MAX;
    	pr->note[i].evidence_level = DBL_MIN;
    	pr->note[i].last_evidence_level = DBL_MIN;
    	pr->note[i].onset_time = DBL_MAX;
    	pr->note[i].state = OFF;
    	pr->note[i].midi_number = i + first_note;

      pr->note[i].oct = (signed char)floor((pr->note[i].midi_number - 60) / 12.0);

    	switch(note_class = pr->note[i].midi_number % 12)
    	{
            case 0:
                strcpy(pr->note[i].name, "C \0");
           		break;
            case 1:
                strcpy(pr->note[i].name, "C#\0");
           		break;
            case 2:
                strcpy(pr->note[i].name, "D \0");
           		break;
            case 3:
                strcpy(pr->note[i].name, "Eb\0");
           		break;
            case 4:
                strcpy(pr->note[i].name, "E \0");
           		break;
            case 5:
                strcpy(pr->note[i].name, "F \0");
           		break;
            case 6:
                strcpy(pr->note[i].name, "F#\0");
           		break;
            case 7:
                strcpy(pr->note[i].name, "G \0");
           		break;
            case 8:
                strcpy(pr->note[i].name, "Ab\0");
           		break;
         	case 9:
            	strcpy(pr->note[i].name, "A \0");
              	break;
            case 10:
                strcpy(pr->note[i].name, "Bb\0");
           		break;
            case 11:
                strcpy(pr->note[i].name, "B \0");
           		break;
    	}
   	}
	LOWEST_PRINTED_NOTE = last_note;
	HIGHEST_PRINTED_NOTE = first_note;
	MIN_PRINTED_NOTE_EVIDENCE_LEVEL = DBL_MAX;
	MAX_PRINTED_NOTE_EVIDENCE_LEVEL = DBL_MIN;

   	return OK;
}


/********************************************************************************************/
/* FREE_PITCH_RANGE                                                                         */
/********************************************************************************************/
void free_pitch_range(Pitch_Range **pr)
{
	unsigned char i;

  free((*pr)->note);
	free(*pr);
}


/********************************************************************************************/
/* TURN_ON_NEW_NOTES                                                                        */
/********************************************************************************************/
/* MAX_DELAY_AFTER_ONSET is assumed to be previously verified */
void turn_on_new_notes(Pitch_Range* pr, const double last_onset_time, const double curr_time)
{
    unsigned char i;

    for (i = 0; i < pr->num_of_notes; i++)
    {
    	if (pr->note[i].state == OFF
    		&& pr->note[i].evidence_level != DBL_MIN
     		&& pr->note[i].last_evidence_time <= curr_time)
   		{
         	turn_on_note(pr, i, last_onset_time);
   		}
    }
}


/********************************************************************************************/
/* TURN_ON_NOTE                                                                             */
/********************************************************************************************/
void turn_on_note(Pitch_Range *pr, const unsigned char note_index,
					const double last_onset_time)
{
	if (pr->note[note_index].state == OFF) {
	 	pr->note[note_index].state = ON;
	 	pr->note[note_index].onset_time = last_onset_time;
	 	pr->current_polyphony++;
	}
}

/********************************************************************************************/
/* TURN_OFF_OBSOLETE_NOTES                                                                  */
/********************************************************************************************/
/* it also prints notes that must be printed before turning them off */
void turn_off_obsolete_notes(Pitch_Range* pr, const double curr_time)
{
	unsigned char i;

    for (i = 0; i < pr->num_of_notes; i++)
        if (pr->note[i].state == ON
        	&& (curr_time - pr->note[i].last_evidence_time) > MAX_GAP_INSIDE_NOTE + TEMPORAL_RESOLUTION)
        {
        	if ( (pr->note[i].last_evidence_time - pr->note[i].onset_time) >= MIN_NOTE_DURATION
        		&& pr->note[i].evidence_level >= MIN_RELATIVE_INTENSITY_LEVEL * MAX_INTENSITY_LEVEL)
        	{
            	print_note(pr->note[i]);
        	}

			turn_off_note(pr, i);
        }
}

/********************************************************************************************/
/* TURN_OFF_NOTE                                                                            */
/********************************************************************************************/
void turn_off_note(Pitch_Range *pr, const unsigned char note_index)
{
	if (pr->note[note_index].state == ON)
	{
		pr->note[note_index].state = OFF;
	    pr->note[note_index].evidence_level = DBL_MIN;
	    pr->note[note_index].last_evidence_level = DBL_MIN;
   	    pr->note[note_index].last_evidence_time = DBL_MAX;
   	    pr->note[note_index].onset_time = DBL_MAX;

	    pr->current_polyphony--;
	}
}


/********************************************************************************************/
/* PRINT_NOTE                                                                               */
/********************************************************************************************/
void print_note(const Note note)
{
	if (PRINT_NOTE_AS == GNUPLOT_ARROWS)
	{
		printf("set arrow from %g, %u, %g to %g, %u, %g arrowstyle 1\n", note.onset_time,
				note.midi_number, note.evidence_level, note.last_evidence_time, note.midi_number, note.evidence_level);
	}
	else if (PRINT_NOTE_AS == MIDIFIABLE_CSV)
	{
    	printf("1, %u, Note_on_c, 0, %u, %u\n",
			    		(unsigned int)((TIME_OFFSET + note.onset_time)*TICS_PER_SEC + 0.5),
			    		note.midi_number,
			    		determine_midi_velocity(note.evidence_level));

    	printf("1, %u, Note_off_c, 0, %u, 0\n",
    			(unsigned int)((TIME_OFFSET + note.last_evidence_time)*TICS_PER_SEC + 0.5),
		    			note.midi_number);

	}
	else if (PRINT_NOTE_AS == HUMAN_READABLE_NOTE_LIST)
	{
    	printf("%.3f\t%s(%2d)\t%.3f\t%.3f\n", note.onset_time, note.name, note.oct,
    											note.last_evidence_time - note.onset_time,
    												note.evidence_level);
	}
	if (note.midi_number < LOWEST_PRINTED_NOTE) {
		LOWEST_PRINTED_NOTE = note.midi_number;
	}
	if (note.midi_number > HIGHEST_PRINTED_NOTE) {
		HIGHEST_PRINTED_NOTE = note.midi_number;
	}
	if (note.evidence_level < MIN_PRINTED_NOTE_EVIDENCE_LEVEL) {
		MIN_PRINTED_NOTE_EVIDENCE_LEVEL = note.evidence_level;
	}
	if (note.evidence_level > MAX_PRINTED_NOTE_EVIDENCE_LEVEL) {
		MAX_PRINTED_NOTE_EVIDENCE_LEVEL = note.evidence_level;
	}
}


/********************************************************************************************/
/* PRINT_PR_STATE                                                                           */
/********************************************************************************************/
void print_pr_state(Pitch_Range* pr, double curr_time)
{
    Note *pr_note_i;
    Note const *pr_last_note = &(pr->note[pr->num_of_notes - 1]);

    for (pr_note_i = &(pr->note[0]); pr_note_i != pr_last_note; pr_note_i++)
        if (pr_note_i->state == ON)
        	printf("%g %u\n", curr_time, pr_note_i->midi_number);
}


/********************************************************************************************/
/* PRINT_PITCH_RANGE                                                                        */
/********************************************************************************************/
void print_pitch_range(Pitch_Range* pr)
{
    unsigned char i;
    unsigned char const num_of_notes = pr->num_of_notes;

    Note *pr_note_i;

    for (i = 0, pr_note_i = &(pr->note[0]);
    	 i != num_of_notes;
      	 i++, pr_note_i++)
  	{
    	printf("%u %u (%g Hz)\n", i, pr->note[i].midi_number, pr->note[i].fund_freq);
    }
}


/********************************************************************************************/
/* PRINT_GNUPLOT_INITIAL_CONFIG                                                             */
/********************************************************************************************/
void print_gnuplot_initial_config(void)
{
	printf("unset arrow\n");

	printf("set ytics 3\n");
	printf("set mytics 3\n");
	printf("set xtics rotate %g, %g\n", INITIAL_PAUSE, 60/TEMPO);
	printf("set grid lw 0.5\n");
	printf("set style arrow 1 nohead lw 2\n");
}


/********************************************************************************************/
/* PRINT_GNUPLOT_FINAL_CONFIG                                                               */
/********************************************************************************************/
void print_gnuplot_final_config(double curr_time)
{
	printf("set xrange [0:%g]\n", curr_time);
	printf("set yrange [%d:%d]\n", LOWEST_PRINTED_NOTE - 1, HIGHEST_PRINTED_NOTE + 1);
	printf("set zrange [%g:%g]\n", MIN_PRINTED_NOTE_EVIDENCE_LEVEL, MAX_PRINTED_NOTE_EVIDENCE_LEVEL);
	printf("splot %g notitle\n", MIN_PRINTED_NOTE_EVIDENCE_LEVEL-1);
}


/********************************************************************************************/
/* PRINT_HUMAN_READABLE_NOTE_LIST_HEADER                                                    */
/********************************************************************************************/
void print_human_readable_note_list_header(void)
{
	printf( "# onset (s)\tduration (s)\n"
			"#\tpitch\t\tintensity\n"
			"\n");
}


/********************************************************************************************/
/* PRINT_CSVMIDI_HEADER                                                                     */
/********************************************************************************************/
void print_csvmidi_header(void)
{
	printf("0, 0, Header, 0, 1, 960\n"
			"1, 0, Start_track\n"
			"1, 0, Time_signature, 4, 2, 24, 8\n");
	printf("1, 0, Tempo, %u\n", TICS_PER_BEAT);
	printf("1, 0, Program_c, 0, %u\n", GM_INSTRUMENT_NUMBER);
	printf("1, 0, Control_c, 0, 7, %u\n", DEFAULT_GM_CHANNEL_VOLUME);
}


/********************************************************************************************/
/* PRINT_CSVMIDI_END_OF_TRACK                                                               */
/********************************************************************************************/
void print_csvmidi_end_of_track(double time)
{
	printf("1, %u, End_track\n", (unsigned int)((TIME_OFFSET + time)*TICS_PER_SEC + 0.5));
}


/********************************************************************************************/
/* PRINT_CSVMIDI_END_OF_FILE                                                                */
/********************************************************************************************/
void print_csvmidi_end_of_file(void)
{
	printf("0, 0, End_of_file\n");
}


/********************************************************************************************/
/* SET_BPM                                                                                  */
/********************************************************************************************/
void set_bpm(const double beats_per_minute)
{
	TICS_PER_BEAT = (unsigned int)(TICS_PER_MIN / beats_per_minute + 0.5);
	TICS_PER_SEC = (double)TICS_PER_MIN / TICS_PER_BEAT * 16;
}


/********************************************************************************************/
/* SET_OFFSET                                                                               */
/********************************************************************************************/
void set_offset(const double offset_in_seconds)
{
	TIME_OFFSET = offset_in_seconds;
}


/********************************************************************************************/
/* SET_INSTRUMENT                                                                           */
/********************************************************************************************/
void set_instrument(const unsigned char instrument_gm_number)
{
	GM_INSTRUMENT_NUMBER = instrument_gm_number;
}


/********************************************************************************************/
/* PRINT_FILE_HEADER                                                                        */
/********************************************************************************************/
void print_file_header(void)
{
	if (PRINT_NOTE_AS == MIDIFIABLE_CSV) {
		print_csvmidi_header();
	} else if (PRINT_NOTE_AS == GNUPLOT_ARROWS) {
		print_gnuplot_initial_config();
	} else if (PRINT_NOTE_AS == HUMAN_READABLE_NOTE_LIST) {
		print_human_readable_note_list_header();
	}
}


/********************************************************************************************/
/* PRINT_FILE_FOOTER                                                                        */
/********************************************************************************************/
void print_file_footer(const double curr_time)
{
	if (PRINT_NOTE_AS == MIDIFIABLE_CSV) {
		print_csvmidi_end_of_track(curr_time);
		print_csvmidi_end_of_file();
	} else if (PRINT_NOTE_AS == GNUPLOT_ARROWS) {
		print_gnuplot_final_config(curr_time);
	}
}


/********************************************************************************************/
/* STIMULATE_NOTE                                                                           */
/********************************************************************************************/
void stimulate_note(Note *note, const double curr_time, const double stimulus_intensity)
{
	/*
	note->last_evidence_time = curr_time;
	if (note->evidence_level < stimulus_intensity) {
		note->evidence_level = stimulus_intensity;
	}
	return;
	*/

	if (stimulus_intensity < MIN_RELATIVE_RESTIMULUS_INTENSITY_LEVEL * MAX_INTENSITY_LEVEL) {
		return;
	}

	if (curr_time - note->onset_time <= MAX_ATTACK_LENGTH) {
		note->last_evidence_time = curr_time;
		note->last_evidence_level = stimulus_intensity;
		if (note->evidence_level < stimulus_intensity) {
			note->evidence_level = stimulus_intensity;
		}
	} else {
		double abs_intens_increase = stimulus_intensity - note->last_evidence_level,
				rel_intens_increase = abs_intens_increase/note->last_evidence_level;

		if (abs_intens_increase > MIN_ATTACK_ABSOLUTE_INCREASE
			&& rel_intens_increase > MIN_ATTACK_RELATIVE_INCREASE
			&& note->last_evidence_time - note->onset_time >= MIN_NOTE_DURATION) { /* termina a nota atual e produz nova nota */
        	if (note->last_evidence_time - note->onset_time >= MIN_NOTE_DURATION
        		&& note->evidence_level >= MIN_RELATIVE_INTENSITY_LEVEL * MAX_INTENSITY_LEVEL)
        	{
            	print_note(*note);
        	}
			note->onset_time = curr_time;
			note->evidence_level = stimulus_intensity;
		}
		note->last_evidence_time = curr_time;
		note->last_evidence_level = stimulus_intensity;
	}
}


/********************************************************************************************/
/* DETERMINE_MIDI_VELOCITY                                                                  */
/********************************************************************************************/
unsigned char determine_midi_velocity(const double evidence_level) {
	if (evidence_level >= MAX_INTENSITY_LEVEL) {
		return 127;
	}
	return 127 * (evidence_level / MAX_INTENSITY_LEVEL);
}


/********************************************************************************************/
/* FREQ2NOTE                                                                                */
/********************************************************************************************/
signed char freq2note(const double freq, const Pitch_Range *pr)
{
	if (freq < pr->note[0].fund_freq || freq > pr->note[pr->num_of_notes-1].fund_freq) {
		return -1;
	}
	return log2(freq / pr->note[0].fund_freq)*12 + 0.5;
}
