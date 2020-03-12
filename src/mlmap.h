#ifndef H_MLMAP
#define H_MLMAP
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "field.h"

void update_V(mpz_t *V, const int num_new, const mpz_t ri);
void evaluate_V(mpz_t rop, const mpz_t *V_in, const int d, const mpz_t *r);

void mlmap_evaluation_N(mpz_t *v, const int n, const uint64 num_terms, const uint64 pivot, const mpz_t* V);

void initialize_beta(mpz_t* betavals, const int d, const mpz_t* z);
void evaluate_beta(mpz_t rop, const mpz_t* z, const mpz_t* r, const int d);

void initialize_alpha(mpz_t *alpha, const int l, const int e);
void evaluate_alpha(mpz_t rop, const uint64 e_in, const mpz_t* r, const int d);

void initialize_tau(mpz_t *tau, const int l, const mpz_t val);
void evaluate_tau(mpz_t rop, const mpz_t val, const mpz_t *r, const int d);

void initialize_tau2(mpz_t *tau, const int l, const mpz_t val);
void evaluate_tau2(mpz_t rop, const mpz_t val, const mpz_t *r, const int d);

void kxi_eval(mpz_t* rop, const uint64 num_terms, const mpz_t* r);

#endif

