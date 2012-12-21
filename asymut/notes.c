#include <math.h>
#include <float.h>
#include <string.h>

#include "notes.h"
#include "aux_fun.h"
#include "numeric_constants.h"
#include "config_constants.h"
#include "config_variables.h"

/********************************************************************************************/
/* INIT_PITCH_RANGE                                                                         */
/********************************************************************************************/
signed char init_pitch_range(Pitch_Range **pitch_range, const unsigned char first_note,
																const unsigned char last_note)
{
    unsigned char i, note_class;
    unsigned char const num_of_notes = last_note - first_note + 1;
    Note *pr_note_i;
    Pitch_Range *pr;
    
    if (last_note < first_note)
    {
        fprintf(stderr, "Error in Pitch_Range initialization: 'first_note' is greater than "
        																	"'last_note'.\n");
    	return ERROR;
   	}
   	
   	if (!*pitch_range) *pitch_range = (Pitch_Range*) emalloc(sizeof(Pitch_Range));
   	
   	pr = *pitch_range;
   	
    pr->note = (Note*) emalloc(num_of_notes * sizeof(Note));
    pr->num_of_notes = num_of_notes;
    pr->current_polyphony = 0;

    for (i = 0, pr_note_i = &(pr->note[0]); i < num_of_notes; i++, pr_note_i++)
    {
        pr_note_i->fund_freq = FREQ_REF_A4 * pow(2, (double)(i + first_note - 69) / 12);
    	pr_note_i->last_evidence_time = (-1) * DBL_MAX;
    	pr_note_i->onset_time = DBL_MAX;
    	pr_note_i->state = OFF;
    	pr_note_i->midi_number = i + first_note;

     	pr_note_i->oct = pr_note_i->midi_number / 12;

    	switch(note_class = pr_note_i->midi_number % 12)
    	{
            case 0:
                strcpy(pr_note_i->name, "C \0");
           		break;
            case 1:
                strcpy(pr_note_i->name, "C#\0");
           		break;
            case 2:
                strcpy(pr_note_i->name, "D \0");
           		break;
            case 3:
                strcpy(pr_note_i->name, "Eb\0");
           		break;
            case 4:
                strcpy(pr_note_i->name, "E \0");
           		break;
            case 5:
                strcpy(pr_note_i->name, "F \0");
           		break;
            case 6:
                strcpy(pr_note_i->name, "F#\0");
           		break;
            case 7:
                strcpy(pr_note_i->name, "G \0");
           		break;
            case 8:
                strcpy(pr_note_i->name, "Ab\0");
           		break;
         	case 9:
            	strcpy(pr_note_i->name, "A \0");
              	break;
            case 10:
                strcpy(pr_note_i->name, "Bb\0");
           		break;
            case 11:
                strcpy(pr_note_i->name, "B \0");
           		break;
    	}
   	}
   	return OK;
}

/********************************************************************************************/
/* TURN_ON_NEW_NOTES                                                                        */
/********************************************************************************************/
/* note that MAX_DELAY_AFTER_ONSET is assumed to be previously verified */
void turn_on_new_notes(Pitch_Range* pr, const double last_onset_time, const double curr_time)
{
    Note *pr_note_i;
    Note const *pr_last_note = &(pr->note[pr->num_of_notes - 1]);
    
    for (pr_note_i = &(pr->note[0]); pr_note_i != pr_last_note; pr_note_i++)
    	if (pr_note_i->state == OFF
     		&& pr_note_i->last_evidence_time == curr_time)
   		{
         	pr_note_i->state = ON;
         	pr_note_i->onset_time = last_onset_time;
         	pr->current_polyphony++;
   		}
}

/********************************************************************************************/
/* TURN_OFF_OBSOLETE_NOTES                                                                  */
/********************************************************************************************/
/* it also prints notes that must be printed before turning them off */
void turn_off_obsolete_notes(Pitch_Range* pr, const double curr_time)
{
    Note *pr_note_i;
    Note const *pr_last_note = &(pr->note[pr->num_of_notes - 1]);

    for (pr_note_i = &(pr->note[0]); pr_note_i != pr_last_note; pr_note_i++)
        if (pr_note_i->state == ON
        	&& (curr_time - pr_note_i->last_evidence_time) > MAX_GAP_INSIDE_NOTE)
        {
        	if ( (pr_note_i->last_evidence_time - pr_note_i->onset_time) >= MIN_NOTE_DURATION)
            	print_note(*pr_note_i);

            pr_note_i->state = OFF;
            pr->current_polyphony--;
        }
}

/********************************************************************************************/
/* PRINT_NOTE                                                                               */
/********************************************************************************************/
void print_note(const Note note)
{
	if (PRINT_NOTE_AS_GNUPLOT_ARROWS)
	{
		printf("set arrow from %g, %u to %g, %u nohead lw 5\n", note.onset_time,
				note.midi_number, note.last_evidence_time, note.midi_number);
	}
	else
	{
    	printf("%lf ; %s(%d) ; %lf\n", note.onset_time, note.name, note.oct,
    											note.last_evidence_time - note.onset_time);
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

