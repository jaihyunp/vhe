#ifndef H_PARAMS
#define H_PARAMS

#include <gmp.h>
#include "typedefs.h"

extern mpz_t PRIME;
extern uint64 logN;
extern uint64 N;
extern mpz_t *ROU;
extern gmp_randstate_t STATE;		// State for random

/* Init */
int init_field(mp_bitcnt_t n, int logN_in);
void free_field();

#endif
