#ifndef _NOTES_H_
#define _NOTES_H_

#include "data_structures.h"

signed char init_pitch_range(Pitch_Range **pr, const unsigned char first_note,
																const unsigned char last_note);
void turn_on_new_notes(Pitch_Range* pr, const double last_onset_time, const double curr_time);
void turn_off_obsolete_notes(Pitch_Range* pr, const double curr_time);
void print_note(const Note note);
void print_pr_state(Pitch_Range* pr, double curr_time);
void print_pitch_range(Pitch_Range* pr);

#endif
