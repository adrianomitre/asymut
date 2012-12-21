#include <stdlib.h>

#include "data_structures.h"
#include "aux_fun.h"

static Postponable_Function *postponed;
static unsigned int num_postponed = 0;

/********************************************************************************************/
/* PREPARE_POSTPONING                                                                       */
/********************************************************************************************/
void prepare_postponing(unsigned int n)
{
	postponed = (Postponable_Function*)emalloc(n*sizeof(Postponable_Function));
}

/********************************************************************************************/
/* FINISH_POSTPONING                                                                        */
/********************************************************************************************/
void finish_postponing(void)
{
	free(postponed);
}

/********************************************************************************************/
/* POSTPONE                                                                                 */
/********************************************************************************************/
void postpone(Postponable_Function fun)
{
	postponed[num_postponed++] = fun;
}

/********************************************************************************************/
/* CALL_POSTPONED                                                                           */
/********************************************************************************************/
void call_postponed(void)
{
	unsigned int i;

	for (i = 0; i < num_postponed; i++)
		postponed[i]();
}
