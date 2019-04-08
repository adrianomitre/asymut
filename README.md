# ASyMuT - Automatic System for Music Transcription (and analysis)

ASyMuT is a system for automatic music transcription (i.e., WAV to MIDI) and analysis.

## Preliminary conventions:

n : natural number  
z : integer number  
r : "real" number (actually decimal <~ floating-point)  
r+ : non-negative "real" number  
s : contiguous string of characters without spaces  
p : one among many predefined strings for the given parameter (multiple-choice or "combo")  
< > : mandatory parameter  
[ ] : optional parameter  
[ ] ... [ ] : variable number of optional parameters  

## Mandatory parameters:

Parameter: `--stdft_win_len <n>`  
Abbreviation: `-w <n>`  
Description: Short-Time Discrete Fourier Transform Window Length (a.k.a Frame Size)  
Associate variable: `unsigned int fft_prop->window->length`  
Observations: It is possible to use the “k”“ multiplier (e.g. “-w 4k” instead of “-w 4096”)  

Parameter: `--input_file <s>`  
Abbreviation: `-i <s>`  
Description: Input File Name  
Associate variable: `char stream->file_name[FILENAME_MAX]`  

## Optional parameters:

Parameter: `--stdft_apod_fun <p> [r] ... [r]`  
Abbreviation: `-a <p>`  
Description: STDFT Apodization Function (a.k.a tapering or smoothing function)  
Associate variable: `Apodization_Function APOD_FUN`  
Possible values: `rectangular`, `triangular`, `hamming`, `hann`, `blackman`, `blackman-harris`, `nuttall3`..`nuttall12`  , `nuttall14`..`nuttall15`, `gaussian`, `hanning-poisson`, `helie_a_w1`, `helie_a_w6`

Parameter: `--stdft_win_step <n>`  
Abbreviation: `-s <n>`  
Description: STDFT Window Step (a.ka. window hop or frame stride)  
Associate variable: `unsigned int fft_prop->step`  

Parameter: `--zero_padding_ratio <n>`  
Abbreviation: `-z <n>`  
Description: Zero Padding Ratio  
Associate variable: `unsigned short fft_prop->zero_padding_ratio`  

Parameter: `--min_abs_f0 <r+>`  
Description: Minimum Absolute F0 (in Hertz)  
Associate variable: `double MIN_ABSOLUTE_F0`  

Parameter: `--min_rel_f0 <r+>`  
Description: Minimum Relative F0 (in Frequency Resolution units)  
Associate variable: `double MIN_F0_TO_FREQ_RES_RATIO`  
Observations: setting to 0 disables this constraint  

Parameter: `--max_abs_f0 <r+>`  
Description: Maximum Absolute F0 (in Hertz)  
Associate variable: `double MAX_ABSOLUTE_F0`  

Parameter: `--max_bands <n>`  
Description: Maximum Critical Bands  
Associate variable: `unsigned char MAX_CRITICAL_BANDS`  

Parameter: `--min_interonset_gap <r+>`  
Description: Minimum Interonset Gap (i.e., Minimum Distance Between Onsets) (in seconds)  
Associate variable: `double MIN_DISTANCE_BETWEEN_ONSETS`  

Parameter: `--min_onset_win <r+>`  
Description: Onset Threshold Window Minimum Length (in seconds)  
Associate variable: `double ONSET_THRESHOLD_ WINDOW_ MIN_LENGTH`  

Parameter: `--max_onset_win <r+>`  
Description: Onset Threshold Window Maximum Length (in seconds)  
Associate variable: `double ONSET_THRESHOLD_ WINDOW_ MAX_LENGTH`  

Parameter: `--percentile <r+>`  
Abbreviation: `-p <r+>`  
Description: Onset Threshold Percentile (in per-unit)  
Associate variable: `double ONSET_THRESHOLD_PERCENTILE`  
Observação:
Rounded to nearest discrete possibility

Parameter: `--scale_factor <r>`  
Abbreviation: `-f <r>`  
Description: Onset Threshold Percentile Scaling Factor  
Associate variable: `double ONSET_THRESHOLD_PERCENTILE_SCALING_FACTOR`  
Observação:
threshold = const_part + scale_factor * percentile(threshold_win)

Parameter: `--const_part <r>`  
Abbreviation: `-c <r>`  
Description: Onset Threshold Constant Part  
Associate variable: `double ONSET_THRESHOLD_CONSTANT_PART`  

Parameter: `--max_delay <r+>`  
Description: Maximum After Onset Delay Before F0 Evidence (in seconds)  
Associate variable: `double MAX_AFTER_ONSET_DELAY_BEFORE_F0_EVIDENCE`  

Parameter: `--max_gap <r+>`  
Description: Maximum F0 Evidence Gap Inside Note (in seconds)  
Associate variable: `double MAX_F0_EVIDENCE_GAP_INSIDE_NOTE`  

Parameter: `--min_duration <r+>`  
Description: Minimum Note Duration (i.e., minimum note evidence time)  
Associate variable: `double MIN_NOTE_DURATION`  

Parameter: `--ref_a4_freq <r+>`  
Description: A4 Reference Frequency (in Hertz)  
Associate variable: ` double FREQ_REF_A4`  

Parameter: `--pow_norm_level`  
Abbreviation: `-n`  
Description: Power Normalization Level (in dB)  
Associate variable: `double POWER_NORMALIZATION_LEVEL`  

Parameter: `--lowest_note <n>`  
Description: Pitch Range Lowest Note  
Associate variable: `unsigned char PITCH_RANGE_LOWEST_NOTE`  

Parameter: `--highest_note <n>`  
Description: Pitch Range Highest Note  
Associate variable: `unsigned char PITCH_RANGE_HIGHEST_NOTE`  

Parameter: `--f0_estimation_method <p> [p] [r] ... [r]`  
Abbreviation: `-m`  
Description: F0 Estimation Method  
Associate variable: `F0_Estimation_Method ESTIMATION_METHOD`  
Possible values: `(max_index, hps, hsc, bw_hsc, fft_fft) [mag|pow|lp|wd], klapuri`  

Parameter: `--unpred_method <p> [r] ... [r]`  
Abbreviation: `-u`  
Description: Unpredictability Estimation Method  
Associate variable: `Unpredictability_Estimation_Method UNPRED_METHOD`  
Possible values: `complex, lp, ilp, sam, ph, new`  

Parameter: `--klap_ic <r+>`  
Description: Klapuri Multiple F0 Estimation Iteration Control parameter  
Associate variable: `double KLAPURI_ITERATION_CONTROL`  

Parameter: `--piano_roll`  
Description: Print Notes as Gnuplot Arrows (piano-roll like)  
Associate variable: `PRINT_NOTES_AS_GNUPLOT_ARROWS`  

## Parameters that change the system purpose
Associate variable: `Asymut_Operation_Mode OPERATION_MODE`  

Parameter: `--transcribe`  
Description: Transcribe (i.e., default system operation mode)  

Parameter: `--print <p> | <p> <p> [p]`  
Description: Prints  
Possible values: `samples, apod_win_time, apod_win_freq, <spec|loc_max> <ph|uph|mag|pow|lp|wd|pn [  doubled_linebreaks]>, unpred,  threshold, bands_response, bands_response_summary, bands_sum, onset, f0_estimate`

Parameter: `--ruminate <p>`  
Description: Ruminates (process but don’t print anything, essencialy for benchmarking)  
Possible values: `spec, unpred, f0_estimate`  

## References

### Undergraduate thesis

#### In Portuguese

[Sistema Automático de Transcrição Melódica](http://www.linux.ime.usp.br/~cef/mac499-05/monografias/amitre/monografia.ps)

### Conference proceedings

#### In Portuguese

[Um Sistema Automático de Transcrição Melódica](http://gsd.ime.usp.br/sbcm/2005/papers/tech-12478.pdf)
conference proceeding

#### In English

##### Slides for the conference proceedings

[A System for the Automatic Transcription of Melodies](https://www.slideshare.net/slideshow/embed_code/key/5eEAomCKJBXb8f)

##### Conference proceedings, focus on the Pitch Determination Algorithm (PDA)

[Accurate and Efficient Fundamental Frequency Determination from Precise Partial Estimates](https://www.ime.usp.br/~mqz/Mitre_AESBR2006.pdf)

##### Related project

PDA extracted (only Hann window + Grandke), no external dependencies, input = samples (not audio files)

[pda](https://github.com/adrianomitre/pda)
