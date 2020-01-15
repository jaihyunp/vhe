#ifndef H_POLY
#define H_POLY

#include "typedefs.h"
#include "gmp.h"

void fft(mpz_t *v, const uint64 n, const mpz_t *c, const uint64 deg);
void ifft(mpz_t *c, const mpz_t *v, const uint64 n);
int polynomial_cmp(const mpz_t *c, const uint64 deg, const mpz_t *v, const uint64 n);
void polynomial_extrapolate(mpz_t val, const mpz_t x, const mpz_t *xs, const mpz_t *ys, const uint64 n);
void polynomial_extrapolate_N(mpz_t val, const mpz_t x, const mpz_t *ys, const uint64 n);

void coeff_add(mpz_t *c, const mpz_t *c1, const mpz_t *c2, const uint64 deg);
void coeff_mult(mpz_t *c, const mpz_t *c1, const mpz_t *c2, const uint64 deg);
void coeff_evaluate(mpz_t val, const mpz_t x, const mpz_t *c, const uint64 deg);

void fourier_add(mpz_t *v, const mpz_t *v1, const mpz_t *v2, const uint64 n);
void fourier_mult(mpz_t *v, const mpz_t *v1, const mpz_t *v2, const uint64 n);
void fourier_extrapolate(mpz_t val, const mpz_t x, const mpz_t *v, const uint64 n);

#endif
