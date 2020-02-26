/*************************************
//Jai Hyun Park
//October 27, 2019.
//Implementation of polynomial rounding circuit using commitments
**************************************/

#ifndef H_VHEAAN
#define H_VHEAAN

//libs
#include <gmp.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "field.h"
#include "mlmap.h"
#include "commit.h"
#include "polynomial.h"

int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL);
//int gkr_cipher_rescale(mpz_t rop, mpz_t *V_d, mpz_t *V_r, mpz_t *V_c, const mpz_t *r_d, const mpz_t *r_r, const mpz_t *r_c, const mpz_t *r_b, const int log_num, const int log_bits, const int bef_bits, const int aft_bits, const mpz_t val);
int gkr_cipher_rescale
(
    const mpz_t *V_r_in, const mpz_t *V_d_in,
    const mpz_t *V_c_in, 
    const mpz_t *r_c, const mpz_t *r_b,
    const int log_num, const int bef_bits, const int aft_bits,
    const mpz_t val
);

int gkr_cipher_mult(
    const mpz_t *V_0r, const mpz_t *V_0c, const mpz_t *V_0d, 
    const mpz_t *V_1l, const mpz_t *V_1r, const mpz_t *V_1c, const mpz_t * V_1d,
    const mpz_t *V_2, 
    const mpz_t *r_0r, const mpz_t *r_0c, const mpz_t *r_0d, const mpz_t *r_0b,
    const mpz_t *r_1l, const mpz_t *r_1r, const mpz_t *r_1c, const mpz_t *r_1d, const mpz_t *r_1b, 
    const mpz_t *r_2, 
    const mpz_t P, const mpz_t EVK1, const mpz_t EVK2, 
    const int log_num, const int bits,
    const mpz_t val
);


//int sum_check_cipher_mult();
//int sum_check_cipher_add();
//int sum_check_cipher_rescale();
//int sum_check_verification();
//
//int sum_check_cipher_mult
//(
//    mpz_t rop, mpz_t *CR1_in, mpz_t *CR2_in, mpz_t *V3, mpz_t *V2, mpz_t *V1l, mpz_t *V1rr, mpz_t *VR1x, mpz_t *VR1y, mpz_t *VR2x, mpz_t *VR2y,
//    const mpz_t ri_in,
//    const mpz_t P, const mpz_t *evk, const mpz_t val,
//    const int logN, const int bit,
//    const mpz_t r, const mpz_t *z,
//    const mpz_t *zR1d, const mpz_t *zR1b, const mpz_t *rR1d, const mpz_t *rR1b,
//    const mpz_t *zR2d, const mpz_t *zR2b, const mpz_t *rR2d, const mpz_t *rR2b,
//    const mpz_t z1, 
//	const mpz_t* commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t* sks, 
//	const mpz_t* commits2, const uint64 lck2, const uint64 key_num2, const uint64 test_comm_num2, const mpz_t* sks2, const mpz_t gen
//);
//
//int sum_check_cipher_cmult(mpz_t rop, mpz_t *V, const mpz_t c, const mpz_t ri, const mpz_t z);
//int sum_check_cipher_add(mpz_t rop, mpz_t *V, const mpz_t ri, const mpz_t r, const mpz_t z);
//
////int sum_check_cipher_rescale(mpz_t rop, mpz_t *C, const mpz_t ri, const int logN, const int l1, const int l2, const mpz_t *rn, const mpz_t *rd, const mpz_t *rb, const mpz_t val, const mpz_t *zn, const mpz_t *zdb);
//
//int sum_check_cipher_rescale(mpz_t rop, mpz_t *C, const mpz_t ri, const int logN, const int l1, const int l2, const mpz_t *rn, const mpz_t *rd, const mpz_t *rb, const mpz_t val, const mpz_t *zn, const mpz_t *zdb,
//	const mpz_t* commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t* sks, const mpz_t gen);
//
//int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL);
//
//
//int sum_check_rounding_fast(mpz_t rop, mpz_t *C, const mpz_t ri, const int d, const int d1, const int d2, const mpz_t *r, const mpz_t* z);
//int sum_check_poly_rounding(mpz_t rop, mpz_t *C, const mpz_t ri, const int ln, const int ld, const int d1, const int d2, const mpz_t *r_n, const mpz_t *r_d, const mpz_t *r_b, const mpz_t val, const mpz_t *z_n, const mpz_t *z_db,
//	const mpz_t* commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t* sks, const mpz_t gen);

//defs
#define ERROR 0x0

#endif
