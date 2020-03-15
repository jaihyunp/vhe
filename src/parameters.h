#ifndef H_PARAMS
#define H_PARAMS

#include <gmp.h>
#include "typedefs.h"
#include "field.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

extern mpz_t PRIME;
extern mpz_t GSIZE;
extern mpz_t GGEN;

extern long WORDS;					// Each element is represented by WORD number of uint64
extern long BITS;					// = WORDS * 64  =  Bit-size of GSIZE
extern long digit_rep;				// Only for printf
extern gmp_randstate_t STATE;		// State for random

extern uint64 logN;
extern uint64 N;

extern mpz_t *ROU;

extern double TIME_PROVER;
extern double TIME_VERIFIER;
extern double TIME_VERIFIER_IOLAYER;
extern double TIME_VERIFIER_GKR;
extern double TIME_PROVER_GKR;
extern double TIME_PROVER_EVAL;

int init_field(mp_bitcnt_t n, int logN_in);
void free_field();

#endif
