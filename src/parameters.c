#include "parameters.h"
#include "field.h"
#include <gmp.h>
#include <time.h>
#include <stdlib.h>

mpz_t PRIME;
uint64 logN, N;
mpz_t *ROU;//2N th roots of unity
gmp_randstate_t STATE;		// State for random


/* Generator */
int is_generator(mpz_t gen, mpz_t p, mpz_t po, const mpz_t max, const mpz_t odd)
{
    int check = 1;
   
    mpz_set(p, odd);
    mod_pow(po, gen, odd);


    while(mpz_cmp(max, p) > 0) {
        
        if(!mpz_cmp_ui(po, 1)) {
            check = 0;
            break;
        }

        mpz_add(p, p, p);
        mod_mult(po, po, po);
    }
    return check;
}

void find_generator(mpz_t gen, const mpz_t odd)
{
    mpz_t p, po, max;
    mpz_inits(p, po, max, NULL);

    mod_sub_ui(max, PRIME, 1);
    mpz_set_ui(gen, 0);
    do {
        mod_add_ui(gen, gen, 1);
    } while(!is_generator(gen, p, po, max, odd));

    mpz_clears(p, po, max, NULL);
}


/* Init */
int init_field(mp_bitcnt_t n, int logN_in)
{
    if((int) n < (logN_in + 1))
        return 0;

    mpz_t odd, gen;
    mpz_inits(PRIME, odd, gen, NULL);

    gmp_randinit_mt(STATE);
    gmp_randseed_ui(STATE, time(NULL));

    printf("====================================================\n");
    printf(" Parameter Settetting\n");
    printf("----------------------------------------------------\n");

    /* N Setting*/
    logN = logN_in;
    N = 1ULL << logN;
    printf("N(logN): %lld(%lld)\n", N, logN);

    /* PRIME Setting */
    do {

        do {
            mpz_rrandomb(odd, STATE, n - logN - 2);
        } while(mpz_probab_prime_p(odd, 50) == 0);

        mpz_mul_2exp(PRIME, odd, logN + 2);
        mpz_add_ui(PRIME, PRIME, 1);

    } while(mpz_probab_prime_p(PRIME, 50) == 0);
    printf("PRIME: %s\n", mpz_get_str(0, 16, PRIME));
   
    /* GEN & ROU Setting */
    find_generator(gen, odd);
    printf("Generator: %s\n", mpz_get_str(0, 16, gen));

    ROU = (mpz_t*) malloc(sizeof(mpz_t) * N * 4);
    for(uint64 i = 0; i < N * 4; i ++)
        mpz_init(ROU[i]);
    mpz_set_ui(ROU[0], 1);
    mod_pow(ROU[1], gen, odd);
    for(uint64 i = 2; i < N * 4; i ++)
        mod_mult(ROU[i], ROU[i - 1], ROU[1]);
    printf("====================================================\n\n\n");

    mpz_clears(odd, gen, NULL);
    return 1;
}

void free_field()
{
    for(uint64 i = 0; i < N * 4; i ++)
        mpz_clear(ROU[i]);
    free(ROU);
    mpz_clears(PRIME, NULL);
}

