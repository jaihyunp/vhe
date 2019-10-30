#include "field.h"
mpz_t PRIME;
void mod(mpz_t rop, const mpz_t x)
{
    mpz_mod(rop, x, PRIME);
}

void mod_add(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_add(rop, x, y);
    mod(rop, rop);
}
void mod_sub(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_sub(rop, x, y);
    mod(rop, rop);
}
void mod_mult(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_mul(rop, x, y);
    mod(rop, rop);
}
void mod_1neg(mpz_t rop, const mpz_t r)
{
    mpz_ui_sub(rop, 1, r);
    mod(rop, rop);
}
void mod_pow(mpz_t rop, const mpz_t x, uint64 b)
{
    uint64 d = 1;
    mpz_t base, res;
    mpz_inits(base, res, NULL);

    mpz_set_ui(res, 1);
    mpz_set(base, x);

    while(d <= b) {
        if(d & b)
            mod_mult(res, res, base);
        d = d << 1;
        if(d <= b)
            mod_mult(base, base, base);
    }

    mod(rop, res);
    mpz_clears(base, res, NULL);
}

void euclidean_algo(const mpz_t u, mpz_t u1, mpz_t u2, mpz_t u3)
{
    mpz_t v1, v2, v3, q, t1, t2, t3, tmp;
    mpz_inits(v1, v2, v3, q, t1, t2, t3, tmp, NULL);

    mpz_set_ui(u1, 1);
    mpz_set_ui(u2, 0);
    mpz_set(u3, u);

    mpz_set_ui(v1, 0);
    mpz_set_ui(v2, 1);
    mpz_set(v3, PRIME);

    do {
        mpz_fdiv_q(q, u3, v3);

        mod_mult(tmp, q, v1);
        mod_sub(t1, u1, tmp);

        mod_mult(tmp, q, v2);
        mod_sub(t2, u2, tmp);

        mod_mult(tmp, q, v3);
        mod_sub(t3, u3, tmp);
        
        mpz_set(u1, v1);
        mpz_set(u2, v2);
        mpz_set(u3, v3);
        
        mpz_set(v1, t1);
        mpz_set(v2, t2);
        mpz_set(v3, t3);
    } while(mpz_sgn(v3) != 0);

    mpz_clears(v1, v2, v3, q, t1, t2, t3, tmp, NULL);
}

void mod_inv(mpz_t rop, const mpz_t a)
{
    mpz_t u1, u2, u3;
    mpz_inits(u1, u2, u3, NULL);

    euclidean_algo(a, u1, u2, u3);

    if(!mpz_cmp_ui(u3, 1))
        mod(rop, u1);
    else
        mpz_set_ui(rop, 0);

    mpz_clears(u1, u2, u3, NULL);
}

void mod_inv_pow2(mpz_t rop, const int k)
{
    mpz_t tmp;
    mpz_init_set_ui(tmp, 2);
    mod_pow(tmp, tmp, k);
    mod_inv(rop, tmp);
    mpz_clear(tmp);
}


