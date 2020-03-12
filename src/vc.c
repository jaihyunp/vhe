#include "vc.h"

int sum_check_verification(mpz_t rop, mpz_t **poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);

    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check\n%s\n%s\n%s\n", str, mpz_get_str(0, digit_rep, tmp), mpz_get_str(0, digit_rep, ri), mpz_get_str(0, digit_rep, PRIME));
        stat = 0;
    }
    
    if(stat) {
        polynomial_extrapolate_N(extrap_val, r[num_rounds - 1], poly[0], degree);
        for(int round = 1; round < num_rounds; round ++) {
            mod_add(tmp, poly[round][0], poly[round][1]);
            if(mpz_cmp(tmp, extrap_val)) {
                printf("Fail sumcheck %s :: %dth round %s %s\n", str, round, mpz_get_str(0, digit_rep, tmp), mpz_get_str(0, digit_rep, extrap_val));
                stat = 0;
                break;
            } else {
                polynomial_extrapolate_N(extrap_val, r[num_rounds - 1 - round], poly[round], degree);
            }
        }
    }

    if(stat)
        mpz_set(rop, extrap_val);
    else
        mpz_set_ui(rop, ERROR);

    mpz_clears(tmp, extrap_val, NULL);

    return stat;
}
