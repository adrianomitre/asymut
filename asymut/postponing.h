#ifndef _POSTPONING_H_
#define _POSTPONING_H_

void prepare_postponing(unsigned int n);
void finish_postponing(void);
void postpone(Postponable_Function fun);
void call_postponed(void);
void test_postponing(void);

#endif /*_POSTPONING_H_ */
