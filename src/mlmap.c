#include "mlmap.h"
void update_V(mpz_t *V, const int num_new, const mpz_t ri)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    for(int i = 0; i < num_new; i ++){
        mod_ui_sub(tmp1, 1, ri);
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
    for(uint64 i = 0; i < n; i ++) {
        mpz_init_set(V[i], V_in[i]);
    }

    for(int round = 0; round < d; round ++) {
        update_V(V, n >> (round + 1), r[d - round - 1]);
    }

    mod(rop, V[0]);
    for(uint64 i = 0; i < n; i ++)
        mpz_clear(V[i]);
    free(V);
}

void mlmap_evaluation_N(mpz_t *v, const int n, const uint64 num_terms, const uint64 pivot, const mpz_t* V)
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

void kxi_eval(mpz_t* rop, const uint64 num_terms, const mpz_t* r) {
	mpz_t tmp;
	mpz_init2(tmp, BITS);
	mod_ui_sub(rop[0], 1, r[0]);
	mpz_set(rop[1], r[0]);
	int k = 1;
	uint64 j;
	if (num_terms > 2) {
		for (uint64 i = 4; i <= num_terms; i <<= 1) {
			j = i;
			do {
				j--;
				if (j & 1)
					mod_mult(rop[j], rop[(j >> 1)], r[k]);
				else {
					mod_ui_sub(tmp, 1, r[k]);
					mod_mult(rop[j], rop[(j >> 1)], tmp);
				}
			} while (j != 0);
			k++;
		}
	}
}

