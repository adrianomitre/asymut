#ifndef _DATA_STRUCTURES_H_
#define _DATA_STRUCTURES_H_

typedef enum {FALSE, TRUE} Boolean;

typedef enum {DEFAULT_OM,
				PRINT_REFINED_F0_ESTIMATE,
				PRINT_GROSS_F0_ESTIMATE,
				PRINT_F0_HARMONIC_SERIES,
				PRINT_VALID_CAND_GROSS_FREQ_AND_AVERAGE_PARTIAL_ENERGY,
				PRINT_VALID_CAND_GROSS_FREQ_AND_TOTAL_ENERGY,
				PRINT_VALID_CAND_ERROR,
				PRINT_INIT_CAND_TOTAL_ENERGY,
				PRINT_INIT_CANDIDATES,
				PRINT_MAX_CAND_FREQ_OM,
				PRINT_MIN_CAND_FREQ_OM,
				PRINT_QUANTILE_OM,
				PRINT_PARAMETERS_OM
				} OperationMode;

typedef struct estimation_parameters_t {
	OperationMode opMode;
	double minAbsFreq;
	double maxAbsFreq;
	double maxRelativeEnergyBelow;  /* Essencially for speeding-up the process */
	double maxRelativeEnergyAbove;
	double minRelativeEnergy;
	double maxRelativeInharmonicity;
} EstimationParameters;

#endif /* _DATA_STRUCTURES_H_ */
