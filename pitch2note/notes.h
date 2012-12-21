#ifndef _NOTES_H_
#define _NOTES_H_

#include "data_structures.h"

signed char init_pitch_range(Pitch_Range **pr, const unsigned char first_note,
																const unsigned char last_note);

void free_pitch_range(Pitch_Range **pr);																
																
void stimulate_note(Note *note, const double curr_time, const double stimulus_intensity);
void turn_on_new_notes(Pitch_Range* pr, const double last_onset_time, const double curr_time);
void turn_on_note(Pitch_Range *pr, const unsigned char note_index,
					const double last_onset_time);
void turn_off_obsolete_notes(Pitch_Range* pr, const double curr_time);
void turn_off_note(Pitch_Range *pr, const unsigned char note_index);

void print_note(const Note note);
void print_pr_state(Pitch_Range* pr, double curr_time);
void print_pitch_range(Pitch_Range* pr);

void print_file_header(void);
void print_file_footer(const double curr_time);

void print_csvmidi_header(void);
void print_csvmidi_end_of_track(double time);
void print_csvmidi_end_of_file(void);

void set_bpm(const double beats_per_minute);
void set_offset(const double offset_in_seconds);
void set_instrument(const unsigned char instrument_gm_number);

unsigned char determine_midi_velocity(const double evidence_level);

signed char freq2note(const double freq, const Pitch_Range *pr);

#endif /* _NOTES_H_ */
