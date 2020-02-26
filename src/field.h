#ifndef H_FIELD
#define H_FIELD

#include "parameters.h"
#include "typedefs.h"
#include "gmp.h"

void mod(mpz_t rop, const mpz_t x);
void mod_set_ui(mpz_t rop, const uint64 x);
void mod_init_set_ui(mpz_t rop, const uint64 x);
void mod_random(mpz_t rop);
void mod_init_random(mpz_t rop);

void mod_add(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_add_ui(mpz_t rop, const mpz_t x, const uint64 y);

void mod_sub(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_sub_ui(mpz_t rop, const mpz_t x, const uint64 y);
void mod_ui_sub(mpz_t rop, const uint64 x, const mpz_t y);
void mod_1neg(mpz_t rop, const mpz_t r); //set rop to be 1-r mod PRIME

void mod_mult(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_mult_ui(mpz_t rop, const mpz_t x, const uint64 y);
void mod_mult_si(mpz_t rop, const mpz_t x, const signed long long y);

void mod_addmul(mpz_t rop, const mpz_t x, const mpz_t y);

void mod_multG(mpz_t rop, const mpz_t x, const mpz_t y);

void mod_pow(mpz_t rop, const mpz_t x, const mpz_t b);
void mod_pow_ui(mpz_t rop, const mpz_t x, const uint64 b);
void mod_ui_pow(mpz_t rop, const uint64 x, const mpz_t b);
void mod_ui_pow_ui(mpz_t rop, const uint64 x, const uint64 b);

void euclidean_algo(const mpz_t u, mpz_t u1, mpz_t u2, mpz_t u3);
void mod_inv(mpz_t rop, const mpz_t a);
void mod_inv_ui(mpz_t rop, const uint64 a);
void mod_inv_pow2(mpz_t rop, const uint64 k);
#endif
