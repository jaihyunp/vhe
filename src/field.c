#include "field.h"
#include <stdio.h>

/*
mpz_t PRIME;
gmp_randstate_t STATE;		// State for random
uint64 logN, N;
mpz_t *ROU;
*/

/* Assign */
void mod(mpz_t rop, const mpz_t x)
{
    mpz_mod(rop, x, PRIME);
}

void mod_set_ui(mpz_t rop, const uint64 x)
{
    char str[20];
    sprintf(str, "%lld", x);
    mpz_set_str(rop, str, 10);
    mod(rop, rop);
}

void mod_init_set_ui(mpz_t rop, const uint64 x)
{
    mpz_init(rop);
    mod_set_ui(rop, x);
}

void mod_random(mpz_t rop)
{
    mpz_urandomm(rop, STATE, PRIME);
}

void mod_init_random(mpz_t rop)
{
    mpz_init(rop);
    mod_random(rop);
}


/* Add */
void mod_add(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_add(rop, x, y);
    mod(rop, rop);
}

void mod_add_ui(mpz_t rop, const mpz_t x, const uint64 y)
{
    mpz_add_ui(rop, x, y);
    mod(rop, rop);
}


/* Sub */
void mod_sub(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_sub(rop, x, y);
    mod(rop, rop);
}

void mod_sub_ui(mpz_t rop, const mpz_t x, const uint64 y)
{
    mpz_sub_ui(rop, x, y);
    mod(rop, rop);
}

void mod_ui_sub(mpz_t rop, const uint64 x, const mpz_t y)
{
    mpz_ui_sub(rop, x, y);
    mod(rop, rop);
}


/* Mult */
void mod_mult(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_mul(rop, x, y);
    mod(rop, rop);
}

void mod_mult_ui(mpz_t rop, const mpz_t x, const uint64 y)
{
    mpz_mul_ui(rop, x, y);
    mod(rop, rop);
}

void mod_mult_si(mpz_t rop, const mpz_t x, const signed long long y)
{
    mpz_mul_si(rop, x, y);
    mod(rop, rop);
}


/* Commit */
//void mod_multG(mpz_t rop, const mpz_t x, const mpz_t y)
//{
//	mpz_mul(rop, x, y);
//	mpz_mod(rop, rop, GSIZE);
//}


/* Pow */
void mod_pow(mpz_t rop, const mpz_t x, const mpz_t b)
{
    mpz_powm(rop, x, b, PRIME);
}

void mod_ui_pow(mpz_t rop, const uint64 x, const mpz_t b)
{
    mpz_t base;
    mod_init_set_ui(base, x);

    mod_pow(rop, base, b);

    mpz_clear(base);
}

void mod_pow_ui(mpz_t rop, const mpz_t x, const uint64 b)
{
    uint64 d = b;
    mpz_t base, res;

    mod_init_set_ui(res, 1);
    mpz_init_set(base, x);

    while(d) {
        if(d & 1)
            mod_mult(res, res, base);
        d >>= 1;
        mod_mult(base, base, base);
    }
    mod(rop, res);
    
    mpz_clears(base, res, NULL);
}

void mod_ui_pow_ui(mpz_t rop, const uint64 x, const uint64 b)
{
    mpz_t base;
    mod_init_set_ui(base, x);

    mod_pow_ui(rop, base, b);

    mpz_clear(base);
}


/* Addmul */ 
void mod_addmul(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_t tmp;
    mpz_init(tmp);
    mod_mult(tmp, x, y);
    mod_add(rop, rop, tmp);
    mpz_clear(tmp);
}

/* Inversion */
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

void mod_inv_ui(mpz_t rop, const uint64 a_in)
{
    mpz_t a;
    mod_init_set_ui(a, a_in);
    mod_inv(rop, a);
    mpz_clear(a);
}

void mod_inv_pow2(mpz_t rop, const uint64 k)
{
    mpz_t tmp;
    mpz_init(tmp);
    mod_ui_pow_ui(tmp, 2, k);
    mod_inv(rop, tmp);
    mpz_clear(tmp);
}


