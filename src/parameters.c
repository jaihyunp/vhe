#include "parameters.h"

mpz_t PRIME;
mpz_t GSIZE;
mpz_t GGEN;
long WORDS, BITS, digit_rep;
uint64 logN, N;
mpz_t *ROU;//2N th roots of unity
gmp_randstate_t STATE;		// State for random
double TIME_PROVER;
double TIME_VERIFIER;

/* Generator */

// Check weather gen is a generator of Z_PRIME or not.
// Note that ord(gen) | 4N*odd, and this script investigates gen^(2^k * odd) for each k with 2^k < 4N.
// Thus, if gen passes this test, its order is 4N*odd, so is a generator of Z_PRIME.
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


// Find a generator of Z_PRIME by exhausitive searching.
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
    mpz_inits(GSIZE, PRIME, odd, gen, NULL);
    digit_rep = 16;

    BITS = n + 1;

    gmp_randinit_mt(STATE);
//    gmp_randseed_ui(STATE, time(NULL));
    gmp_randseed_ui(STATE, 0);

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
            do {
                mpz_rrandomb(odd, STATE, n - logN - 2);
            } while (mpz_probab_prime_p(odd, 50) == 0); // set odd to be prime
            
            mpz_mul_2exp(PRIME, odd, logN + 2);
            mpz_add_ui(PRIME, PRIME, 1);

        } while (mpz_probab_prime_p(PRIME, 50) == 0); // set PRIME = 4N*odd + 1 to be prime

		mpz_add(GSIZE, PRIME, PRIME);
		mpz_add_ui(GSIZE, GSIZE, 1);
	} while (mpz_probab_prime_p(GSIZE, 50) == 0); // set GSIZE = 2*PRIME + 1 to be prime

    printf("PRIME: %s\n", mpz_get_str(0, digit_rep, PRIME));
    printf("GSIZE: %s\n", mpz_get_str(0, digit_rep, GSIZE));
   
    /* GEN & ROU Setting */
    find_generator(gen, odd);
    printf("Gen of Z_PRIME: %s\n", mpz_get_str(0, digit_rep, gen));

    ROU = (mpz_t*) malloc(sizeof(mpz_t) * N * 4);
    for (uint64 i = 0; i < N * 4; i ++)
        mpz_init(ROU[i]);
    mpz_set_ui(ROU[0], 1);
    mod_pow(ROU[1], gen, odd);
    for (uint64 i = 2; i < N * 4; i ++)
        mod_mult(ROU[i], ROU[i - 1], ROU[1]);


    /* Commit Setting */
	mpz_t genTest;
	mpz_init2(genTest, BITS);
	mpz_init2(GGEN, BITS);

	do {
		do {
			mpz_urandomb(GGEN, STATE, BITS);
			mpz_mod(GGEN, GGEN, GSIZE);
			mod_multG(genTest, GGEN, GGEN);
		} while (mpz_cmp_ui(genTest, 1) == 0);			// repeat untill order of gen is not 2

		mpz_powm(genTest, GGEN, PRIME, GSIZE);

	} while (mpz_cmp_ui(genTest, 1) != 0);				// repeat untill order of gen is PRIME

	printf("gen of Z_GSIZE: %s\n", mpz_get_str(NULL, digit_rep, GGEN));

    printf("====================================================\n");

    mpz_clears(gen, odd, NULL);
    TIME_PROVER = 0;
    TIME_VERIFIER = 0;
    return 1;
}

void free_field()
{
    for(uint64 i = 0; i < N * 4; i ++)
        mpz_clear(ROU[i]);
    free(ROU);
    mpz_clears(PRIME, GSIZE, GGEN, NULL);
}

