/*************************************
//Jai Hyun Park
//October 27, 2019.
//Implementation of polynomial rounding circuit using commitments
**************************************/

#include "mlmap.h"

void extrapolate(mpz_t rop, const mpz_t *vec, const int n, const mpz_t r)
{
    mpz_t mult, tmp, inv;
    mpz_inits(mult, tmp, inv, NULL);

    mpz_set_ui(rop, 0);
    for(int i = 0; i < n; i ++) {
        mpz_set_ui(mult, 1);
        for(int j = 0; j < n; j ++) {
            if(i != j){
                mpz_set_si(tmp, i - j);
                mod_inv(inv, tmp);
                mpz_sub_ui(tmp, r, (unsigned int) j);
                mod_mult(mult, mult, tmp);
                mod_mult(mult, mult, inv);
            }
        }
        mod_mult(tmp, mult, vec[i]);
        mod_add(rop, rop, tmp);
    }

    mpz_clears(mult, tmp, inv, NULL);
}


void update_V(mpz_t *V, const int num_new, const mpz_t ri)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    for(int i = 0; i < num_new; i ++){
        mod_1neg(tmp1, ri);
        mod_mult(tmp1, tmp1, V[i]);

        mod_mult(tmp2, ri, V[i + num_new]);

        mod_add(V[i], tmp1, tmp2);
    }
    mpz_clears(tmp1, tmp2, NULL);
}

void evaluate_V(mpz_t rop, const mpz_t *V_in, const int d, const mpz_t *r)
{
    uint64 n = 1ULL << d;
    mpz_t *V = (mpz_t*) malloc(sizeof(mpz_t) << d);
    for(uint64 i = 0; i < n; i ++)
        mpz_init_set(V[i], V_in[i]);

    for(int round = 0; round < d; round ++) {
        update_V(V, n >> (round + 1), r[d - round - 1]);
    }

    mod(rop, V[0]);
    for(uint64 i = 0; i < n; i ++)
        mpz_clear(V[i]);
    free(V);
}

void linear_evaluation(mpz_t *v, const int n, const uint64 num_terms, const int pivot, const mpz_t* V)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    mpz_set(v[0], V[pivot]);
    mpz_set(v[1], V[pivot + (num_terms >> 1)]);

    for(int i = 2; i < n; i ++) {
        mpz_mul_si(tmp1, v[1], i);
        mpz_mul_si(tmp2, v[0], 1 - i);
        mod_add(v[i], tmp1, tmp2);
    }

    mpz_clears(tmp1, tmp2, NULL);
}


