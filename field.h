#ifndef H_FIELD
#define H_FIELD
#include "vheaan.h"
extern mpz_t PRIME;
void mod(mpz_t rop, const mpz_t x);
void mod_add(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_sub(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_mult(mpz_t rop, const mpz_t x, const mpz_t y);
void mod_1neg(mpz_t rop, const mpz_t r);
void mod_pow(mpz_t rop, const mpz_t x, uint64 b);
void euclidean_algo(const mpz_t u, mpz_t u1, mpz_t u2, mpz_t u3);
void mod_inv(mpz_t rop, const mpz_t a);
void mod_inv_pow2(mpz_t rop, const int k);
#endif
