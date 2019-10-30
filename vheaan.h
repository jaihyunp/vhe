/*************************************
//Jai Hyun Park
//October 27, 2019.
//Implementation of polynomial rounding circuit using commitments
**************************************/

#ifndef H_VHEAAN
#define H_VHEAAN

//typedefs
typedef unsigned long long uint64;

//libs
#include <gmp.h>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "field.h"
#include "mlmap.h"

int sum_check_cipher_cmult(mpz_t rop, mpz_t *V, const mpz_t c, const mpz_t ri, const mpz_t z);
int sum_check_cipher_add(mpz_t rop, mpz_t *V, const mpz_t ri, const mpz_t r, const mpz_t z);
int sum_check_cipher_rescale(mpz_t rop, mpz_t *C, const mpz_t ri, const int logN, const int l1, const int l2, const mpz_t *rn, const mpz_t *rd, const mpz_t *rb, const mpz_t val, const mpz_t *zn, const mpz_t *zdb);


int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL);


int sum_check_poly_rounding(mpz_t rop, mpz_t *C, const mpz_t ri, const int ln, const int ld, const int d1, const int d2, const mpz_t *r_n, const mpz_t *r_d, const mpz_t *r_b, const mpz_t val, const mpz_t *z_n, const mpz_t *z_db);
int sum_check_rounding_fast(mpz_t rop, mpz_t *C, const mpz_t ri, const int d, const int d1, const int d2, const mpz_t *r, const mpz_t* z);

//defs
#define ERROR 0x0

#endif
