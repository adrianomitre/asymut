#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

typedef enum {FALSE, TRUE} Boolean;
typedef enum {OFF, ON} State;
typedef enum {DOWN, UP} Flag;
typedef enum {HUMAN_READABLE_NOTE_LIST, GNUPLOT_ARROWS, MIDIFIABLE_CSV} Note_Printing_Style;

typedef struct t_note
{
	char name[3];
    signed char oct;
    unsigned char midi_number;

    double fund_freq;

    State state;

    double onset_time; /* in seconds */
    double last_evidence_time; /* in seconds */

	double evidence_level; /* actually should be refactored to max_evidence_level */
	double last_evidence_level;

/*	double last_salience_level; */ /* even if not sounding, note has a small level of presence
 * 									  (something like its probability of being present) */
}
Note;

typedef struct t_pitch_range
{
    unsigned char num_of_notes;
    unsigned char current_polyphony;
  
    Note* note;
}
Pitch_Range;

#endif /* _DATA_STRUCTURES_H_ */
