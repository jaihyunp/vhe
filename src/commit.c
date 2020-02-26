#include "commit.h"


// Given gen of Group G, generate 'key_num' of keys ('pks' for Prover, 'sks' for Verifier) for commitment
void commit_keygen(mpz_t* pks, mpz_t* sks, const uint64 key_num) {
	for (uint64 i = 0; i < key_num; ++i) {
		mpz_init2(sks[i], BITS);
		mpz_init2(pks[i], BITS);
		mpz_urandomb(sks[i], STATE, BITS);
		mod(sks[i], sks[i]);
		mpz_powm(pks[i], GGEN, sks[i], GSIZE);

		//printf("SK: %s\n", mpz_get_str(NULL, digit_rep, sks[i]));
		//printf("PK: %s\n", mpz_get_str(NULL, digit_rep, pks[i]));
	}
}


// Given 'input' vector over F_P of size 'input_row_num' * 'key_num', and 'pks' for commit, output 'commits'. 
// The 'input_row_num' * 2**'log_bit_size' is the total number of output commits. 
// Each input is commited as 2**'log_bit_size' binaries. at 2**'log_bit_size'-th bit to from LSB to MSB Left to Right
void commit_commit(mpz_t* commits, const mpz_t* input, const uint64 key_num, const uint64 input_row_num, const long log_bit_size, const mpz_t* pks) {
	mpz_t idx;
	mpz_t test;
	mpz_init2(test, BITS);
	uint64 bit_size = (1 << log_bit_size);
	
	for (uint64 i = 0; i < input_row_num; ++i) {
		mpz_init2(idx, BITS);
		mpz_add_ui(idx, idx, 1);						
		// mpz_mul_2exp(idx, idx, bit_size-1);				// set idx = (1 << bit_size)

		for (uint64 k = 0; k < bit_size; ++k) {
			mpz_init2(commits[bit_size*i + k], BITS);
			mpz_add_ui(commits[bit_size*i + k], commits[bit_size*i + k], 1);		// set commits[bit_size*i+k] = 1

			for (uint64 j = 0; j < key_num; ++j) {
				mpz_and(test, input[i + j*input_row_num], idx);
				if (mpz_sgn(test) == 1)
					mod_multG(commits[bit_size*i + k], commits[bit_size*i + k], pks[j]);
			}

			mpz_mul_2exp(idx, idx, 1);
			// mpz_divexact_ui(idx, idx, 2);
		}

	}
}


// Given 'evalpt', vector 'output' of size 'key_num' 
// which is a matrixt-vector mult of BD matrix 'input' of size 'key_num' by 'input_row_num' 
// with vector 'evalpt' of size 'input_row_num'*2**'log_bit_size' 
// Vf checks if 'commits' dots 'evalpt' equals 'sks' dots 'output'
void commit_open(mpz_t* output, const mpz_t* input, const mpz_t* evalpts, const mpz_t* commits, const uint64 key_num, const uint64 input_row_num, const long log_bit_size, const mpz_t* sks) {

	clock_t st = clock();

// Prover calculates the result in plain state.
	mpz_t idx;
	mpz_t test;
	mpz_init2(test, BITS);
	uint64 bit_size = (1 << log_bit_size);


	for (uint64 j = 0; j < key_num; ++j) {
		mpz_init2(output[j], BITS);

		for (uint64 i = 0; i < input_row_num; ++i) {
			mpz_init2(idx, BITS);
			mpz_add_ui(idx, idx, 1);
			// mpz_mul_2exp(idx, idx, bit_size - 1);				// set idx = (1 << bit_size)

			for (uint64 k = 0; k < bit_size; ++k) {
				mpz_and(test, input[i + j*input_row_num], idx);
				if (mpz_sgn(test) == 1)
					mod_add(output[j], output[j], evalpts[bit_size*i + k]);

				mpz_mul_2exp(idx, idx, 1);
				//mpz_divexact_ui(idx, idx, 2);
			}
		}
	}

	std::cout << "Prover time is: " << (double)((double)clock() - st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;

	st = clock();

// Verifier checks the validity of output
	mpz_t tmp;
	mpz_init2(tmp, BITS);
	mpz_init2(test, BITS);
	mpz_add_ui(test, test, 1);
	for (uint64 i = 0; i < input_row_num*bit_size; ++i) {
		mpz_powm(tmp, commits[i], evalpts[i], GSIZE);
		mod_multG(test, test, tmp);
	}

	mpz_t test2;
	mpz_init2(test2, BITS);
	for (uint64 i = 0; i < key_num; ++i) {
		mod_mult(tmp, output[i], sks[i]);
		mod_add(test2, test2, tmp);
	}
	mpz_powm(tmp, GGEN, test2, GSIZE);
	
	if (mpz_cmp(test, tmp) == 0)
		std::cout << "Correct commit!" << std::endl;
	else
		std::cout << "Wrong commit!" << std::endl;

	std::cout << "Verifier time is: " << (double)((double)clock() - st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
}



// Given 'input' binary vector of size 'input_num' * 'key_num', and 'pks' for commit, output 'commits'. 
void commit_commit_binary(mpz_t* commits, const mpz_t* input, const uint64 key_num, const uint64 input_num, const mpz_t* pks) {
	mpz_t test, idx;
	mpz_init2(test, BITS);
	mpz_init2(idx, BITS);
	mpz_add_ui(idx, idx, 1);

	for (uint64 i = 0; i < input_num; ++i) {
		//printf("%d \n", i);
		mpz_init2(commits[i], BITS);
		mpz_add_ui(commits[i], commits[i], 1);		// set commits[bit_size*i+k] = 1
		for (uint64 j = 0; j < key_num; ++j) {
			//printf("%d \n", j);
			mpz_and(test, input[i + j*input_num], idx);
			if (mpz_sgn(test) == 1)
				mod_multG(commits[i], commits[i], pks[j]);
		}
	}
}


// Given 'evalpt', vector 'output' of size 'key_num' 
// which is a matrixt-vector mult of Binary matrix 'input' of size 'key_num' by 'input_num' 
// with vector 'evalpt' of size 'input_num'
// Vf checks if 'commits' dots 'evalpt' equals 'sks' dots 'output'
void commit_open_binary(mpz_t* output, const mpz_t* input, const mpz_t* evalpts, const mpz_t* commits, const uint64 key_num, const uint64 input_num, const mpz_t* sks) {

	clock_t st = clock();

	// Prover calculates the result in plain state.
	mpz_t idx, test;
	mpz_init2(test, BITS);
	mpz_init2(idx, BITS);
	mpz_add_ui(idx, idx, 1);
	
	for (uint64 j = 0; j < key_num; ++j) {
		mpz_init2(output[j], BITS);

		for (uint64 i = 0; i < input_num; ++i) {
			mpz_and(test, input[i + j*input_num], idx);
			if (mpz_sgn(test) == 1)
				mod_add(output[j], output[j], evalpts[i]);
		}
	}

	std::cout << "Prover time is: " << (double)((double)clock() - st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;

	st = clock();

	// Verifier checks the validity of output
	mpz_t tmp;
	mpz_init2(tmp, BITS);
	mpz_init2(test, BITS);
	mpz_add_ui(test, test, 1);
	for (uint64 i = 0; i < input_num; ++i) {
		mpz_powm(tmp, commits[i], evalpts[i], GSIZE);
		mod_multG(test, test, tmp);
	}

	mpz_t test2;
	mpz_init2(test2, BITS);
	for (uint64 i = 0; i < key_num; ++i) {
		mod_mult(tmp, output[i], sks[i]);
		mod_add(test2, test2, tmp);
	}
	mpz_powm(tmp, GGEN, test2, GSIZE);

		
	if (mpz_cmp(test, tmp) == 0)
		std::cout << "Correct commit!" << std::endl;
	else
		std::cout << "Wrong commit!" << std::endl;

	std::cout << "Verifier time is: " << (double)((double)clock() - st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
}


void commit_test() {

	clock_t test_st;
	gmp_randinit_mt(STATE);
	gmp_randseed_ui(STATE, time(NULL));

	//init_finte_gp(gen, 1);			// Initialize finite group for Commitment, and set the generator 'gen' of Group.

	uint64 test_key_num = (1 << 11);
	uint64 test_in_num = (1 << 14);
	long test_log_bits = 8;
	uint64 test_comm_num = (1 << test_log_bits) * test_in_num / test_key_num;

	mpz_t* test_input = (mpz_t*)malloc((test_in_num) * sizeof(mpz_t));
	mpz_t* test_sk = (mpz_t*)malloc((test_key_num) * sizeof(mpz_t));
	mpz_t* test_pk = (mpz_t*)malloc((test_key_num) * sizeof(mpz_t));
	mpz_t* test_outs = (mpz_t*)malloc((test_key_num) * sizeof(mpz_t));
	mpz_t* commits = (mpz_t*)malloc((test_comm_num) * sizeof(mpz_t));
	mpz_t* test_evalpts = (mpz_t*)malloc((test_comm_num) * sizeof(mpz_t));

	test_st = clock();
	commit_keygen(test_pk, test_sk, test_key_num);					// Key Generation for Commitment
	std::cout << "KeyGen time is: " << (double)((double)clock() - test_st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;

	for (uint64 i = 0; i < test_in_num; ++i) {							// Generate Random Inputs for Test
		mpz_init2(test_input[i], BITS);
		mpz_urandomb(test_input[i], STATE, BITS);
		mod(test_input[i], test_input[i]);
	}

	for (uint64 i = 0; i < test_comm_num; ++i) {							// Generate Random Evalpts for Test
		mpz_init2(test_evalpts[i], BITS);
		mpz_urandomb(test_evalpts[i], STATE, BITS);
		mod(test_evalpts[i], test_evalpts[i]);
	}

	test_st = clock();
	commit_commit(commits, test_input, test_key_num, test_in_num / test_key_num, test_log_bits, test_pk);			// Prover Commits! 
	std::cout << "Commit time is: " << (double)((double)clock() - test_st) / (CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;

	std::cout << "start open" << std::endl;

	commit_open(test_outs, test_input, test_evalpts, commits, test_key_num, test_in_num / test_key_num, test_log_bits, test_sk); // Prover Evaluate Results and Decommits! 



}