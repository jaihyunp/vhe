#ifndef H_COMMIT
#define H_COMMIT
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "field.h"
#include "parameters.h"

void commit_keygen(mpz_t* pks, mpz_t* sks, const uint64 key_num);

void commit_commit(mpz_t* commits, const mpz_t* input, const uint64 key_num, const uint64 input_row_num, const mpz_t* pks);

void commit_open(mpz_t* output, const mpz_t* input, const mpz_t* evalpts, const mpz_t* commits, const uint64 key_num, const uint64 input_row_num, const mpz_t* sks);

void commit_commit_binary(mpz_t* commits, const mpz_t* input, const uint64 key_num, const uint64 input_num, const mpz_t* pks);

void commit_open_binary(mpz_t* output, const mpz_t* input, const mpz_t* evalpts, const mpz_t* commits, const uint64 key_num, const uint64 input_num, const mpz_t* sks);


//void commit_test();

#endif
