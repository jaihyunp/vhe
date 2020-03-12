#ifndef H_VC
#define H_VC

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "field.h"
#include "mlmap.h"
#include "commit.h"
#include "polynomial.h"

int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL);

#endif
