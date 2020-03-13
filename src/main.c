#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include "mlmap.h"
#include "vc.h"
#include "commit.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int main_FX_mult(int argc, char **argv);
int main_FXoverPhi_mult(int argc, char **argv);
int main_RXoverPhi_mult(int argc, char **argv);

int main(int argc, char **argv)
{
    return main_FX_mult(argc, argv);
//    return main_FXoverPhi_mult(argc, argv);
//    return main_RXoverPhi_mult(argc, argv);
}


int main_RXoverPhi_mult(int argc, char **argv)
{
    if(argc != 7) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num] [bits] [log_ck0_bits] [log_ck1].\n");
        return 0;
    }
    clock_t st;
    double elapsed_time;
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), num = 1 << log_num, stat = 1,
    		log_bits = 0, bits = atoi(argv[4]), log_ck0 = atoi(argv[5]), log_ck1 = atoi(argv[6]), ubits;
    while ((bits - 1) >> ++ log_bits);
    ubits = 1 << log_bits;

    uint64 CK0num = 1 << log_ck0;											// number of commit keys0
    uint64 CM0num = 1 << (log_num + logN_in + 2 + log_bits - log_ck0);		// number of commits0
    uint64 CK1num = 1 << log_ck1;											// number of commit keys1
    uint64 CM1num = 1 << (log_num + logN_in + 1 - log_ck1);					// number of commits1


    init_field(bits_of_prime, logN_in); //init field

    //// variables for commit0
    mpz_t *sks0 = (mpz_t *) malloc(sizeof(mpz_t) * CK0num),
    		*pks0 = (mpz_t *) malloc(sizeof(mpz_t) * CK0num),
			*CK0_out = (mpz_t *) malloc(sizeof(mpz_t) * CK0num),
			*CK0_in = (mpz_t *) malloc(sizeof(mpz_t) * 4 * N * num * ubits),
			*tmp_C0r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2 + logN + log_bits)),
			*tmp_rC0r = (mpz_t *) malloc(sizeof(mpz_t) * log_ck0),
			*commits0 = (mpz_t *) malloc(sizeof(mpz_t) * CM0num),
			*C0r_evalpts = (mpz_t *) malloc(sizeof(mpz_t) * CM0num);
    for (uint64 i = 0; i < CK0num; i++)
    	mpz_inits(sks0[i], pks0[i], CK0_out[i], 0);
    for (uint64 i = 0; i < 4 * N * num * ubits; i++)
    	mpz_init(CK0_in[i]);
    for (uint64 i = 0; i < log_num + 2 + logN + log_bits; i++)
    	mpz_init(tmp_C0r[i]);
    for (int i = 0; i < log_ck0; i++)
    	mpz_init(tmp_rC0r[i]);
    for (uint64 i = 0; i < CM0num; i++)
        mpz_inits(commits0[i], C0r_evalpts[i], 0);

    //// variables for commit1
    mpz_t *sks1 = (mpz_t *) malloc(sizeof(mpz_t) * CK1num),
    		*pks1 = (mpz_t *) malloc(sizeof(mpz_t) * CK1num),
			*CK1_out = (mpz_t *) malloc(sizeof(mpz_t) * CK1num),
			*CK1_in = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N * num),
			*tmp_C1r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN)),
			*tmp_rC1r = (mpz_t *) malloc(sizeof(mpz_t) * log_ck1),
			*commits1 = (mpz_t *) malloc(sizeof(mpz_t) * CM1num),
			*C1r_evalpts = (mpz_t *) malloc(sizeof(mpz_t) * CM1num);
    for (uint64 i = 0; i < CK1num; i++)
		mpz_inits(sks1[i], pks1[i], CK1_out[i], 0);
    for (uint64 i = 0; i < 2 * N * num; i++)
    	mpz_init(CK1_in[i]);
    for (uint64 i = 0; i < log_num + 1 + logN; i++)
    	mpz_init(tmp_C1r[i]);
    for (int i = 0; i < log_ck1; i++)
      	mpz_init(tmp_rC1r[i]);
    for (uint64 i = 0; i < CM1num; i++)
    	mpz_inits(commits1[i], C1r_evalpts[i], 0);

    mpz_t val, V0z, V1z, V2z, V3r, C1r, C0r, tmp, rir, b1r, b2r, b3r, tr, t2r, a1qr, a4qr, COMMIT_R_1, COMMIT_R_0, BITMAX, Q, SIGN;
    mpz_inits(val, V0z, V1z, V2z, V3r, C1r, C0r, tmp, rir, b1r, b2r, b3r, tr, t2r, a1qr, a4qr, COMMIT_R_1, COMMIT_R_0, BITMAX, Q, SIGN, 0);
    mpz_t *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp3 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(tmp1[i], tmp2[i], tmp3[i], 0);

    mpz_t *beta1 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *beta2 = (mpz_t *) malloc(sizeof(mpz_t) * N),
          *beta3 = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
          *alpha_4Q = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
          *alpha_1Q = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
          *tau_2N = (mpz_t *) malloc(sizeof(mpz_t) * N * 2),
          *tau_N = (mpz_t *) malloc(sizeof(mpz_t) * N);
    for (int i = 0; i < num; i ++)
        mpz_init(beta1[i]);
    for (int i = 0; i < (int) N; i ++)
        mpz_init(beta2[i]);
    for (int i = 0; i < ubits * 4; i ++)
        mpz_inits(beta3[i], alpha_4Q[i], alpha_1Q[i], 0);
    for (int i = 0; i < (int) N * 2; i ++)
        mpz_inits(tau_2N[i], 0);
    for (int i = 0; i < (int) N; i ++)
        mpz_inits(tau_N[i], 0);
    mpz_t **poly_V3_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * log_num),
          **poly_C1_V1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1)),
          **poly_C1_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1)),
          **poly_C0_V0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
          **poly_C0_V1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
          **poly_C0_bits = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2));
    for (int i = 0; i < log_num + (int) logN + log_bits + 2; i ++) {
        poly_C0_V0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_C0_V1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_C0_bits[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);

        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_C0_V0[i][j], 0);
            mpz_init_set_ui(poly_C0_V1[i][j], 0);
            mpz_init_set_ui(poly_C0_bits[i][j], 0);
        }
    }
    for (int i = 0; i < log_num; i ++) {
        poly_V3_V2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly_V3_V2[i][j], 0);
    }
    for (int i = 0; i < log_num + (int) logN + 1; i ++) {
        poly_C1_V1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_C1_V2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_C1_V1[i][j], 0);
            mpz_init_set_ui(poly_C1_V2[i][j], 0);
        }
    }
    mpz_t b1[4], b2[4], b3[4], v30[4], v31[4], c1[4], c0[4], t[4], t2[4], a4q[4], a1q[4], sign[4];
    for (int i = 0; i < 4; i ++)
        mpz_inits(b1[i], b2[i], b3[i], v30[i], v31[i], c1[i], c0[i], t[i], t2[i], a4q[i], a1q[i], sign[i], 0);

    mpz_t **v_3 = (mpz_t **) malloc(sizeof(mpz_t *) * num * 2),
          **v_2 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_1 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_0 = (mpz_t **) malloc(sizeof(mpz_t *) * num),

          *C_1 = (mpz_t *) malloc(sizeof(mpz_t) * num * N * 2),
          *C_0 = (mpz_t *) malloc(sizeof(mpz_t) * num * N * ubits * 4),

          *V_3 = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *V_2 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_1 = (mpz_t *) malloc(sizeof(mpz_t) * num * N),
          *V_0 = (mpz_t *) malloc(sizeof(mpz_t) * num);

    mod_ui_pow_ui(Q, 2, bits);
    for (int i = 0 ; i < num; i ++) {
        mpz_inits(V_2[i], V_0[i], 0);
        for (uint64 j = 0; j < N; j ++)
            mpz_init(V_1[i * N + j]);
        v_2[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        v_1[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        v_0[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N * 2; j ++)
            mpz_init(v_2[i][j]);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_3[i]);
        v_3[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N; j ++) {
            mpz_init(v_3[i][j]);
            mpz_urandomm(v_3[i][j], STATE, Q);
        }
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_3 be the array of polynomials of degree N. (1)
    }
    for (int i = 0; i < num * (int) N * 2; i ++)
        mpz_init(C_1[i]);
    for (int i = 0; i < num * (int) N * ubits * 4; i ++)
        mpz_init(C_0[i]);
    mod_ui_pow_ui(BITMAX, 2, 4 * ubits - 1);
    printf("1.1. Memory allocated and Input generated\n");

    // Vf generates the commit keys (1-1)
    st = clock();
    commit_keygen(pks0, sks0, CK0num);
    commit_keygen(pks1, sks1, CK1num);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("1.2. Verifier genertaed commit key.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);



    // Pv evaluates the circuit and commit the required values. (2)
    // Pv evaluates the circuit over polynomial rings. (2-1)
    // Assume that Vf gives and takes the **coefficients** of input polynomials; hence, Pv should do (i)fft.
    st = clock();
    for (int i = 0; i < num; i ++) {
        fft(tmp1, 2 * N, v_3[2 * i], N);
        fft(tmp2, 2 * N, v_3[2 * i + 1], N);
        fourier_mult(tmp3, tmp1, tmp2, 2 * N);
        ifft(v_2[i], tmp3, 2 * N);
        for (uint64 j = 0; j < 2 * N; j ++)
            mpz_set(C_1[i * 2 * N + j], v_2[i][j]);
        for (uint64 j = 0; j < N; j ++) {
            mod_sub(v_1[i][j], v_2[i][j], v_2[i][j + N]);
            mod_add(v_1[i][j], v_1[i][j], BITMAX);  // Is this safe?
            for (int k = 0; k < ubits * 4; k ++) { //additional cost?
                if (mpz_tstbit(v_1[i][j], k))
                    mpz_set_ui(C_0[(i * N + j) * ubits * 4 + k], 1);
                else
                    mpz_set_ui(C_0[(i * N + j) * ubits * 4 + k], 0);
            }
            mpz_mod(v_0[i][j], v_1[i][j], Q);
        }
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("2.1. Circuit evaluated.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // Pv commits the values (2-2).
    st = clock();
    for (int i = 0; i < num; i ++)
    	for (uint64 j = 0; j < N; j ++)
    		for (int k = 0; k < ubits * 4; k ++)
    			mpz_set(CK0_in[(i * N + j) * ubits * 4 + k], C_0[(i * N + j) * ubits * 4 + k]);

    commit_commit_binary(commits0, CK0_in, CK0num, CM0num, pks0);
    printf("Prover commits0 \n");

    for (int i = 0; i < num; i ++)
    		for (uint64 j = 0; j < 2 * N; j ++)
    		    mpz_set(CK1_in[i * 2 * N + j], C_1[i * 2 * N + j]);

    commit_commit(commits1, CK1_in, CK1num, CM1num, pks1);
    printf("Prover commits1 \n");
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("2.2. Prover committed the values.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);



    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3)
    // Vf also should do the same procedure but only on the input and output layer.
    st = clock();
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_0[i], val, v_0[i], N); // Both Pv and Vf
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N); // Both Pv and Vf
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    TIME_VERIFIER += elapsed_time;
    printf("3.1.1. Construct OPR circuit on input and output layers.\n  Both +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);
    st = clock();
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N); // Pv
    for (int i = 0; i < num; i ++)
        for (uint64 j = 0; j < N; j ++)
            mpz_set(V_1[i * N + j], v_1[i][j]); // Pv
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("3.1.2. Construct OPR circuit of the entire circuit.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // GKR Protocol (4)
    // Vf fixes the randomness, and precompute the desired values of input and output layer. (4-1)
    mpz_t *z1 = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *z2 = (mpz_t *) malloc(sizeof(mpz_t) * logN),
          *z3 = (mpz_t *) malloc(sizeof(mpz_t) * (log_bits + 2)),
          *r1 = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *r2 = (mpz_t *) malloc(sizeof(mpz_t) * (logN + 1)),
          *r3 = (mpz_t *) malloc(sizeof(mpz_t) * (log_bits + 2)),
          *tmp_r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + 1 + log_bits + 2));
    for (int i = 0; i < log_num; i ++)
        mod_init_random(z1[i]);
    for (int i = 0; i < log_num + 1; i ++)
        mod_init_random(r1[i]);
    for (int i = 0; i < (int) logN; i ++)
        mod_init_random(z2[i]);
    for (int i = 0; i < (int) logN + 1; i ++)
        mod_init_random(r2[i]);
    for (int i = 0; i < log_bits + 2; i ++) {
        mod_init_random(z3[i]);
        mod_init_random(r3[i]);
    }
    for (int i = 0; i < log_num + 1 + (int) logN + 1 + log_bits + 2; i ++)
        mpz_init(tmp_r[i]);
    printf("4.1. Randomness fixed.\n");


    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r1[i + 1]);

    ////////// Commit1
    ////////// save random point for commited value evaluation.
    for (int i = 0; i < (int) logN + 1 + log_num; ++i)
        mpz_set(tmp_C1r[i], tmp_r[logN + log_num - i]);				// tmp_Cr => reverse order
	for (int i = 0; i < log_ck1; ++i)
		mpz_set(tmp_rC1r[i], tmp_C1r[log_ck1 - 1 - i]);  				// tmp_rCr => order


    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(tmp_r[i], r3[i]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(tmp_r[i + log_bits + 2], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + log_bits + 2], r1[i + 1]);
//    evaluate_V(COMMIT_R_0, C_0, log_num + logN + log_bits + 2, tmp_r);


    ////////// Commit0
    ////////// save random point for commited value evaluation.
    for (int i = 0; i < (int) logN + 2 + log_bits + log_num; ++i)
        mpz_set(tmp_C0r[i], tmp_r[logN + 1 + log_bits + log_num - i]);			// tmp_Cr => reverse order

	for (int i = 0; i < log_ck0; ++i)
		mpz_set(tmp_rC0r[i], tmp_C0r[log_ck0 - 1 - i]);  				// tmp_rCr => order

	////////// Commit done




    // Vf evaluates the MLE of input and output layers at the chosen point. (4-2)
    st = clock();
    evaluate_V(V0z, V_0, log_num, z1);
    evaluate_V(V3r, V_3, log_num + 1, r1);
    printf("....Prover Precomputed.\n");
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.2. Verifier precompute the correct answer.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Pv generates the proof for the circuit. (4-3)
    // Pv gives poly_V3_V2[][][], V_3[0], V_3[1], and C1r to Vf.
    st = clock();
    initialize_beta(beta1, log_num, z1);
    initialize_beta(beta2, logN, z2);
    initialize_beta(beta3, log_bits + 2, z3);
    initialize_alpha(alpha_4Q, log_bits + 2, ubits * 4);
    initialize_alpha(alpha_1Q, log_bits + 2, bits);
    initialize_tau(tau_2N, logN + 1, val);
    initialize_tau(tau_N, logN, val);

    // 0 ~ log_num-1 th round.
    uint64 num_terms = 1ULL << log_num;
    for (int round = 0; round < log_num; round ++) {
        for (uint64 p = 0; p < num_terms / 2; p ++) {
            mlmap_evaluation_N(b1, 4, num_terms, p, beta1);
            mlmap_evaluation_N(v30, 4, num_terms << 1, p << 1, V_3);
            mlmap_evaluation_N(v31, 4, num_terms << 1, (p << 1) + 1, V_3);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b1[i], v30[i]);
                mod_mult(tmp, tmp, v31[i]);
                mod_add(poly_V3_V2[round][i], poly_V3_V2[round][i], tmp);
            }

            for (uint64 q = 0; q < N * 2; q ++) {
                mlmap_evaluation_N(c1, 4, num_terms << (logN + 1), (p << (logN + 1)) + q, C_1);
                mlmap_evaluation_N(t2, 4, 0, q, tau_2N);

                for (int i = 0; i < 4; i ++) {
                    mod_mult(tmp, b1[i], c1[i]);
                    mod_mult(tmp, tmp, t2[i]);
                    mod_add(poly_C1_V2[round][i], poly_C1_V2[round][i], tmp);
                }

                if (q < N) {
                    mlmap_evaluation_N(t, 4, 0, q, tau_N);
                    mlmap_evaluation_N(b2, 4, 0, q, beta2);

                    for (int i = 0; i < 4; i ++) {
                        mod_mult(tmp, b1[i], b2[i]);
                        mod_mult(tmp, tmp, c1[i]);
                        mod_add(poly_C1_V1[round][i], poly_C1_V1[round][i], tmp);
                    }

                    for (int e = 0; e < ubits * 4; e ++) {
                        mlmap_evaluation_N(a4q, 4, 0, e, alpha_4Q);
                        mlmap_evaluation_N(a1q, 4, 0, e, alpha_1Q);
                        mlmap_evaluation_N(b3, 4, 0, e, beta3);
                        mlmap_evaluation_N(c0, 4, num_terms << (logN + log_bits + 2), (((p << logN) + q) << (log_bits + 2)) + e, C_0);

                        for (int i = 0; i < 4; i ++) {
                            mod_mult(tmp, b1[i], b2[i]);
                            mod_mult(tmp, tmp, a4q[i]);
                            mod_mult(tmp, tmp, c0[i]);
                            mod_add(poly_C0_V1[round][i], poly_C0_V1[round][i], tmp);

                            mod_mult(tmp, b1[i], t[i]);
                            mod_mult(tmp, tmp, a1q[i]);
                            mod_mult(tmp, tmp, c0[i]);
                            mod_add(poly_C0_V0[round][i], poly_C0_V0[round][i], tmp);

                            mod_sub_ui(tmp, c0[i], 1);
                            mod_mult(tmp, tmp, c0[i]);
                            mod_mult(tmp, tmp, b1[i]);
                            mod_mult(tmp, tmp, b2[i]);
                            mod_mult(tmp, tmp, b3[i]);
                            mod_add(poly_C0_bits[round][i], poly_C0_bits[round][i], tmp);
                        }
                    }
                } else {
                    mlmap_evaluation_N(t, 4, 0, q - N, tau_N);
                    mlmap_evaluation_N(b2, 4, 0, q - N, beta2);

                    for (int i = 0; i < 4; i ++) {
                        mod_sub(tmp, PRIME, c1[i]);
                        mod_mult(tmp, tmp, b1[i]);
                        mod_mult(tmp, tmp, b2[i]);
                        mod_add(poly_C1_V1[round][i], poly_C1_V1[round][i], tmp);
                    }
                }
            }
        }
        num_terms >>= 1;
        update_V(V_3, num_terms << 1, r1[log_num + 1 - round - 1]);
        update_V(C_1, num_terms << (logN + 1), r1[log_num + 1 - round - 1]);
        update_V(C_0, num_terms << (logN + log_bits + 2), r1[log_num + 1 - round - 1]);
        update_V(beta1, num_terms, r1[log_num + 1 - round - 1]);
    }
    mlmap_evaluation_N(b1, 4, num_terms, 0, beta1);

    // log_num th round.
    for (int i = 0; i < 4; i ++) {
        mpz_set_si(sign[i], 1 - 2 * i);
        mod(sign[i], sign[i]);
    }
    num_terms = 1ULL << (logN + 1);
    for (uint64 q = 0; q < num_terms / 2; q ++) {
        mlmap_evaluation_N(c1, 4, num_terms, q, C_1);
        mlmap_evaluation_N(t2, 4, num_terms, q, tau_2N);
        mlmap_evaluation_N(b2, 4, 0, q, beta2);
        mlmap_evaluation_N(t, 4, 0, q, tau_N);

        for (int i = 0; i < 4; i ++) {
            mod_mult(tmp, b1[i], c1[i]);
            mod_mult(tmp, tmp, t2[i]);
            mod_add(poly_C1_V2[log_num][i], poly_C1_V2[log_num][i], tmp);

            mod_mult(tmp, sign[i], c1[i]);
            mod_mult(tmp, tmp, b1[i]);
            mod_mult(tmp, tmp, b2[i]);
            mod_add(poly_C1_V1[log_num][i], poly_C1_V1[log_num][i], tmp);
        }
    }
    num_terms >>= 1;
    update_V(C_1, num_terms, r2[logN]);
    update_V(tau_2N, num_terms, r2[logN]);
    for (int i = 0; i < 4; i ++) {
        mod_ui_sub(sign[i], 1, r2[logN]);
        mod_sub(sign[i], sign[i], r2[logN]);
    }

    // log_num+1 ~ log_num+logN (or log_num ~ log_num+logN-1)th round.
    for (int round = 0; round < (int) logN; round ++) {
        for (uint64 q = 0; q < num_terms / 2; q ++) {
            mlmap_evaluation_N(c1, 4, num_terms, q, C_1);
            mlmap_evaluation_N(b2, 4, num_terms, q, beta2);
            mlmap_evaluation_N(t, 4, num_terms, q, tau_N);
            mlmap_evaluation_N(t2, 4, num_terms, q, tau_2N);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b1[i], c1[i]);
                mod_mult(tmp, tmp, t2[i]);
                mod_add(poly_C1_V2[round + log_num + 1][i], poly_C1_V2[round + log_num + 1][i], tmp);

                mod_mult(tmp, sign[i], c1[i]);
                mod_mult(tmp, tmp, b1[i]);
                mod_mult(tmp, tmp, b2[i]);
                mod_add(poly_C1_V1[round + log_num + 1][i], poly_C1_V1[round + log_num + 1][i], tmp);
            }

            for (int e = 0; e < ubits * 4; e ++) {
                mlmap_evaluation_N(a4q, 4, 0, e, alpha_4Q);
                mlmap_evaluation_N(a1q, 4, 0, e, alpha_1Q);
                mlmap_evaluation_N(b3, 4, 0, e, beta3);
                mlmap_evaluation_N(c0, 4, num_terms << (log_bits + 2), (q << (log_bits + 2)) + e, C_0);

                for (int i = 0; i < 4; i ++) {
                    mod_mult(tmp, b1[i], b2[i]);
                    mod_mult(tmp, tmp, a4q[i]);
                    mod_mult(tmp, tmp, c0[i]);
                    mod_add(poly_C0_V1[round + log_num][i], poly_C0_V1[round + log_num][i], tmp);

                    mod_mult(tmp, b1[i], t[i]);
                    mod_mult(tmp, tmp, a1q[i]);
                    mod_mult(tmp, tmp, c0[i]);
                    mod_add(poly_C0_V0[round + log_num][i], poly_C0_V0[round + log_num][i], tmp);

                    mod_sub_ui(tmp, c0[i], 1);
                    mod_mult(tmp, tmp, c0[i]);
                    mod_mult(tmp, tmp, b1[i]);
                    mod_mult(tmp, tmp, b2[i]);
                    mod_mult(tmp, tmp, b3[i]);
                    mod_add(poly_C0_bits[round + log_num][i], poly_C0_bits[round + log_num][i], tmp);
                }
            }
        }
        num_terms >>= 1;
        update_V(C_1, num_terms, r2[logN - round - 1]);
        update_V(C_0, num_terms << (log_bits + 2), r2[logN - round - 1]);
        update_V(beta2, num_terms, r2[logN - round - 1]);
        update_V(tau_N, num_terms, r2[logN - round - 1]);
        update_V(tau_2N, num_terms, r2[logN - round - 1]);
    }
    mlmap_evaluation_N(b2, 4, num_terms, 0, beta2);
    mlmap_evaluation_N(t2, 4, num_terms, 0, tau_2N);
    mlmap_evaluation_N(t, 4, num_terms, 0, tau_N);
    mlmap_evaluation_N(c1, 4, num_terms, 0, C_1);

    // last layers.
    num_terms = 1ULL << (log_bits + 2);
    for (int round = 0; round < log_bits + 2; round ++) {
        for (uint64 e = 0; e < num_terms / 2; e ++) {
            mlmap_evaluation_N(a4q, 4, num_terms, e, alpha_4Q);
            mlmap_evaluation_N(a1q, 4, num_terms, e, alpha_1Q);
            mlmap_evaluation_N(b3, 4, num_terms, e, beta3);
            mlmap_evaluation_N(c0, 4, num_terms, e, C_0);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b1[i], b2[i]);
                mod_mult(tmp, tmp, a4q[i]);
                mod_mult(tmp, tmp, c0[i]);
                mod_add(poly_C0_V1[round + log_num + logN][i], poly_C0_V1[round + log_num + logN][i], tmp);

                mod_mult(tmp, b1[i], t[i]);
                mod_mult(tmp, tmp, a1q[i]);
                mod_mult(tmp, tmp, c0[i]);
                mod_add(poly_C0_V0[round + log_num + logN][i], poly_C0_V0[round + log_num + logN][i], tmp);

                mod_sub_ui(tmp, c0[i], 1);
                mod_mult(tmp, tmp, c0[i]);
                mod_mult(tmp, tmp, b1[i]);
                mod_mult(tmp, tmp, b2[i]);
                mod_mult(tmp, tmp, b3[i]);
                mod_add(poly_C0_bits[round + log_num + logN][i], poly_C0_bits[round + log_num + logN][i], tmp);
            }
        }
        num_terms >>= 1;
        update_V(C_0, num_terms, r3[log_bits + 2 - round - 1]);
        update_V(beta3, num_terms, r3[log_bits + 2 - round - 1]);
        update_V(alpha_4Q, num_terms, r3[log_bits + 2 - round - 1]);
        update_V(alpha_1Q, num_terms, r3[log_bits + 2 - round - 1]);
    }
    mpz_set(C0r, C_0[0]);
    mpz_set(C1r, C_1[0]);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("4.3. Proof generated.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // Vf verifies the proof. (4-4)
    st = clock();
    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(tmp_r[i], r3[i]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(tmp_r[i + log_bits + 2], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + log_bits + 2], r1[i + 1]);
    mpz_set_ui(tmp, 0);
    stat = stat ? sum_check_verification(rir, poly_C0_bits, tmp, 4, log_num + logN + log_bits + 2, tmp_r, "GKR C0->bits") : 0;
    if (stat) {
        for (int i = 0; i < log_num; i ++)
            mpz_set(tmp_r[i], r1[i + 1]);
        evaluate_beta(b1r, z1, tmp_r, log_num);
        evaluate_beta(b2r, z2, r2, logN);
        evaluate_beta(b3r, z3, r3, log_bits + 2);
        mod_sub_ui(tmp, C0r, 1);
        mod_mult(tmp, tmp, C0r);
        mod_mult(tmp, tmp, b1r);
        mod_mult(tmp, tmp, b2r);
        mod_mult(tmp, tmp, b3r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR C0->bits Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        } else {
            printf("....GKR C0->bits sumcheck verified.\n");
        }
    }

    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(tmp_r[i], r3[i]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(tmp_r[i + log_bits + 2], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + log_bits + 2], r1[i + 1]);
    stat = stat ? sum_check_verification(rir, poly_C0_V0, V0z, 4, log_num + logN + log_bits + 2, tmp_r, "GKR V0->C0") : 0;
    if (stat) {
        evaluate_alpha(a1qr, bits, r3, log_bits + 2);
        evaluate_alpha(a4qr, 4 * ubits, r3, log_bits + 2);
        evaluate_tau(tr, val, r2, logN);
        evaluate_tau(t2r, val, r2, logN + 1);
        mod_mult(tmp, b1r, tr);
        mod_mult(tmp, tmp, a1qr);
        mod_mult(tmp, tmp, C0r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V0->C0 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        } else {
            printf("....GKR V0->C0 sumcheck verified.\n");
        }
    }

    mod_add(V1z, poly_C0_V1[0][0], poly_C0_V1[0][1]);
    stat = stat ? sum_check_verification(rir, poly_C0_V1, V1z, 4, log_num + logN + log_bits + 2, tmp_r, "GKR V1->C0") : 0;
    if (stat) {
        mod_mult(tmp, b1r, b2r);
        mod_mult(tmp, tmp, a4qr);
        mod_mult(tmp, tmp, C0r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V1->C0 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        } else {
            printf("....GKR V1->C0 sumcheck verified.\n");
        }
    }

    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r1[i + 1]);
    mod_sub(V1z, V1z, BITMAX);
    stat = stat ? sum_check_verification(rir, poly_C1_V1, V1z, 4, log_num + logN + 1, tmp_r, "GKR V1->C1") : 0;
    if (stat) {
        mod_mult(tmp, b1r, b2r);
        mod_mult(tmp, tmp, sign[0]);
        mod_mult(tmp, tmp, C1r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V1->C1 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        } else {
            printf("....GKR V1->C1 sumcheck verified.\n");
        }
    }

    mod_add(V2z, poly_C1_V2[0][0], poly_C1_V2[0][1]);
    stat = stat ? sum_check_verification(rir, poly_C1_V2, V2z, 4, log_num + logN + 1, tmp_r, "GKR V2->C1") : 0;
    if (stat) {
        mod_mult(tmp, b1r, t2r);
        mod_mult(tmp, tmp, C1r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V2->C1 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        } else {
            printf("....GKR V2->C1 sumcheck verified.\n");
        }
    }


    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i], r1[i + 1]);
    stat = stat ? sum_check_verification(rir, poly_V3_V2, V2z, 4, log_num, tmp_r, "GKR V2->V3") : 0;
    if (stat) {
        printf("....GKR V2->V3 sumcheck verified.\n");
        mod_mult(tmp, b1r, V_3[0]);
        mod_mult(tmp, tmp, V_3[1]);
        if (mpz_cmp(tmp, rir)) {
            printf("GKR V2->V3 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
               mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        }
    }

    if (stat) {
        mod_1neg(rir, r1[0]);
        mod_mult(rir, rir, V_3[0]);
        mod_mult(tmp, r1[0], V_3[1]);
        mod_add(rir, rir, tmp); // rir = (1-r0)V3(R0)+r0V3(R1) (= V3(Rr))
        printf("..GKR all verified.\n..Returns polys, C0r (%s), C1r (%s) and rir (%s).\n",
                mpz_get_str(0, digit_rep, C0r),
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, rir));
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.4. Proof verified.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);



    // Vf checks the consistencies: rir and V3r / C1r, C0r and commit (5)

    // Vf checks whether C1r is consistent with the commit. (5-1)
    st = clock();
    kxi_eval(C0r_evalpts, CM0num, &(tmp_C0r[log_ck0]));
    kxi_eval(C1r_evalpts, CM1num, &(tmp_C1r[log_ck1]));
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.1.1. commit: kxi eval.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Open commit
    commit_open_binary(CK0_out, CK0_in, C0r_evalpts, commits0, CK0num, CM0num, sks0);
    printf("Opening Commit0 done!\n");
    commit_open(CK1_out, CK1_in, C1r_evalpts, commits1, CK1num, CM1num, sks1);

    st = clock();
	evaluate_V(COMMIT_R_0, CK0_out, log_ck0, tmp_rC0r);
	evaluate_V(COMMIT_R_1, CK1_out, log_ck1, tmp_rC1r);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.1.3. commit: evaluate_V.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);
    ////////// Commit

    st = clock();
    if(mpz_cmp(COMMIT_R_0, C0r)) {
        printf("Commit Fail: Inconsistent with the committed 0 value: %s %s\n",
                mpz_get_str(0, digit_rep, C0r),
                mpz_get_str(0, digit_rep, COMMIT_R_0));
        stat = 0;
    } else {
        printf("....Consistent with the committed 0 value.\n");
    }
    if(mpz_cmp(COMMIT_R_1, C1r)) {
        printf("Commit Fail: Inconsistent with the committed 1 value: %s %s\n",
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, COMMIT_R_1));
        stat = 0;
    } else {
        printf("....Consistent with the committed 1 value.\n");
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.1.4. commit: check consistency\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // Vf checks whether rir is consistent with the precomputed V3r. (5-2).
    st = clock();
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    } else {
        printf("....Consistent with the original input.\n");
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.2. Check consistency with the input\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    if (stat)
        printf("..Consistency verified.\n");

    printf("Result: %s\n", stat ? "PASS" : "FAIL");
    printf("Total time:\nVf: %f ms.\nPv: %f ms.\n", TIME_VERIFIER, TIME_PROVER);
    return stat;
}


//Implementation of Verifiable Multiplication over Z_p[X]/(X^N + 1), where p is prime and N is power of two.
int main_FXoverPhi_mult(int argc, char **argv)
{
    if(argc != 5) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num] [logCKnum].\n");
        return 0;
    }
    clock_t st;
    double elapsed_time;
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), log_ck = atoi(argv[4]), num = 1 << log_num, stat = 1;
    uint64 CKnum = 1ULL << log_ck;							// number of commit keys
    uint64 CMnum = 1ULL << (log_num + logN_in + 1 - log_ck);		// number of commits
    init_field(bits_of_prime, logN_in); //init field


    //// variables for commit
    mpz_t *sks = (mpz_t *) malloc(sizeof(mpz_t) * CKnum),
    			*pks = (mpz_t *) malloc(sizeof(mpz_t) * CKnum),
			*CK_out = (mpz_t *) malloc(sizeof(mpz_t) * CKnum),
			*CK_in = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N * num),
			*tmp_Cr = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN)),
			*tmp_rCr = (mpz_t *) malloc(sizeof(mpz_t) * log_ck),
			*commits = (mpz_t *) malloc(sizeof(mpz_t) * CMnum),
			*Cr_evalpts = (mpz_t *) malloc(sizeof(mpz_t) * CMnum);
    for (uint64 i = 0; i < CKnum; i++)
    		mpz_inits(sks[i], pks[i], CK_out[i], 0);
    for (uint64 i = 0; i < 2 * N * num; i++)
    		mpz_init(CK_in[i]);
    for (uint64 i = 0; i < log_num + 1 + logN; i++)
        	mpz_init(tmp_Cr[i]);
    for (int i = 0; i < log_ck; i++)
            	mpz_init(tmp_rCr[i]);
    for (uint64 i = 0; i < CMnum; i++)
        	mpz_inits(commits[i], Cr_evalpts[i], 0);


    mpz_t val, V1z, V2z, V3r, C1r, tmp, rir, br, tr, t2r, COMMIT_R;
    mpz_inits(val, V1z, V2z, V3r, C1r, tmp, rir, br, tr, t2r, COMMIT_R, 0);
    mpz_t *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp3 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(tmp1[i], tmp2[i], tmp3[i], 0);

    mpz_t *beta = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *tau_2N = (mpz_t *) malloc(sizeof(mpz_t) * N * 2),
          *tau_N = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
    for (int i = 0; i < num; i ++)
        mpz_init(beta[i]);
    for (int i = 0; i < (int) N * 2; i ++)
        mpz_inits(tau_2N[i], tau_N[i], 0);
    mpz_t **poly_V3_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * log_num),
          **poly_C1_V1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1)),
          **poly_C1_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1));
    for (int i = 0; i < log_num; i ++) {
        poly_V3_V2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly_V3_V2[i][j], 0);
    }
    for (int i = 0; i < log_num + (int) logN + 1; i ++) {
        poly_C1_V1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_C1_V2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_C1_V1[i][j], 0);
            mpz_init_set_ui(poly_C1_V2[i][j], 0);
        }
    }

    mpz_t b[4], v30[4], v31[4], c[4], t[4], t2[4];
    for (int i = 0; i < 4; i ++)
        mpz_inits(b[i], v30[i], v31[i], c[i], t[i], t2[i], 0);


    mpz_t **v_3 = (mpz_t **) malloc(sizeof(mpz_t *) * num * 2),
          **v_2 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_1 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          *C_1 = (mpz_t *) malloc(sizeof(mpz_t) * num * N * 2),
          *V_3 = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *V_2 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_1 = (mpz_t *) malloc(sizeof(mpz_t) * num);

    for (int i = 0 ; i < num; i ++) {
        mpz_inits(V_2[i], V_1[i], 0);
        v_2[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        v_1[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N * 2; j ++)
            mpz_inits(v_2[i][j], v_1[i][j], 0);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_3[i]);
        v_3[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_3[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_3 be the array of polynomials of degree N. (1)
    }
    for (int i = 0; i < num * (int) N * 2; i ++)
        mpz_init(C_1[i]);

    printf("1.1. Memory allocated and Input Generated.\n");

    ////////// Vf commitkey gen (1-1)
    st = clock();
    commit_keygen(pks, sks, CKnum);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("1.2. Verifier genertaed commit key.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Pv evaluates the circuit and commit the required values. (2)
    // Pv evaluates the circuit over polynomial rings. (2-1)
    // Assume that Vf gives and takes the **coefficients** of input polynomials; hence, Pv should do (i)fft.
    st = clock();
    for (int i = 0; i < num; i ++) {
        fft(tmp1, 2 * N, v_3[2 * i], N);
        fft(tmp2, 2 * N, v_3[2 * i + 1], N);
        fourier_mult(tmp3, tmp1, tmp2, 2 * N);
        ifft(v_2[i], tmp3, 2 * N);
        for (uint64 j = 0; j < 2 * N; j ++)
            mpz_set(C_1[i * 2 * N + j], v_2[i][j]);
        for (uint64 j = 0; j < N; j ++)
            mod_sub(v_1[i][j], v_2[i][j], v_2[i][j + N]);
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("2.1. Prover evaluates the circuit.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    ////////// Pv commits the values (2-2)
    st = clock();
    for (int i = 0; i < num; i ++)
    		for (uint64 j = 0; j < 2 * N; j ++)
    		    mpz_set(CK_in[i * 2 * N + j], v_2[i][j]);

    commit_commit(commits, CK_in, CKnum, CMnum, pks);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("2.2. Prover committed the values.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);



    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3)
    // Vf also should do the same procedure but only on the input and output layer.
    mod_random(val);
    st = clock();
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_1[i], val, v_1[i], N); // Both Pv and Vf
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N); // Both Pv and Vf
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    TIME_VERIFIER += elapsed_time;
    printf("3.1.1. Construct one point reduced in-out layers.\n  Both +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    st = clock();
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N); // Pv
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("3.1.2. Construct one point reduced circuit.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // GKR Protocol (4)
    // Vf fixes the randomness, and precompute the desired values of input and output layer. (4-1)
    mpz_t *z = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + 1)),
          *tmp_r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + 1));
    for (int i = 0; i < log_num; i ++)
        mod_init_random(z[i]);
    for (int i = 0; i < log_num + 1 + (int) logN + 1; i ++) {
        mod_init_random(r[i]);
        mpz_init(tmp_r[i]);
    }
    printf("4.1. Randomness fixed.\n"); 


    ////////// Commit
    ////////// save random point for commited value evaluation.
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_Cr[logN + log_num - i], r[i]);					// tmp_Cr => reverse order
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_Cr[log_num - i - 1], r[i + 1 + logN + 1]);

	for (int i = 0; i < log_ck; ++i)
		mpz_set(tmp_rCr[i], tmp_Cr[log_ck - 1 - i]);  				// tmp_rCr => order



    // Vf evaluates the MLE of input and output layers at the chosen point. (4-2)
    st = clock();
    evaluate_V(V1z, V_1, log_num, z);
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(tmp_r[i], r[i + logN + 1]);
    evaluate_V(V3r, V_3, log_num + 1, tmp_r);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.2. Verifier computed the correct answers.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Pv generates the proof for the circuit. (4-3) 
    // Pv gives poly_V3_V2[][][], V_3[0], V_3[1], and C1r to Vf.
    st = clock();
    initialize_beta(beta, log_num, z);
    initialize_tau2(tau_N, logN + 1, val);
    initialize_tau(tau_2N, logN + 1, val);
    uint64 num_terms = 1ULL << log_num;
    for (int round = 0; round < log_num; round ++) {
        for (uint64 p = 0; p < num_terms / 2; p ++) {
            mlmap_evaluation_N(b, 4, num_terms, p, beta);
            mlmap_evaluation_N(v30, 4, num_terms << 1, p << 1, V_3);
            mlmap_evaluation_N(v31, 4, num_terms << 1, (p << 1) + 1, V_3);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b[i], v30[i]);
                mod_mult(tmp, tmp, v31[i]);
                mod_add(poly_V3_V2[round][i], poly_V3_V2[round][i], tmp);
            }

            for (uint64 q = 0; q < N * 2; q ++) {
                mlmap_evaluation_N(c, 4, num_terms << (logN + 1), (p << (logN + 1)) + q, C_1);
                mlmap_evaluation_N(t, 4, 0, q, tau_N);
                mlmap_evaluation_N(t2, 4, 0, q, tau_2N);

                for (int i = 0; i < 4; i ++) {
                    mod_mult(rir, b[i], c[i]);
                    mod_mult(tmp, rir, t2[i]);
                    mod_add(poly_C1_V2[round][i], poly_C1_V2[round][i], tmp);
                    
                    mod_mult(tmp, rir, t[i]);
                    mod_add(poly_C1_V1[round][i], poly_C1_V1[round][i], tmp);
                }
            }
        }
        num_terms >>= 1;
        update_V(V_3, num_terms << 1, r[log_num + 1 + logN + 1 - round - 1]);
        update_V(C_1, num_terms << (logN + 1), r[log_num + 1 + logN + 1 - round - 1]);
        update_V(beta, num_terms, r[log_num + 1 + logN + 1 - round - 1]);
    }

    mlmap_evaluation_N(b, 4, num_terms, 0, beta);
    num_terms = 1ULL << (logN + 1);
    for (int round = log_num; round < log_num + (int) logN + 1; round ++) {
        for (uint64 q = 0; q < num_terms / 2; q ++) {
            mlmap_evaluation_N(c, 4, num_terms, q, C_1);
            mlmap_evaluation_N(t, 4, num_terms, q, tau_N);
            mlmap_evaluation_N(t2, 4, num_terms, q, tau_2N);

            for (int i = 0; i < 4; i ++) {
                mod_mult(rir, b[i], c[i]);
                mod_mult(tmp, rir, t2[i]);
                mod_add(poly_C1_V2[round][i], poly_C1_V2[round][i], tmp);
                
                mod_mult(tmp, rir, t[i]);
                mod_add(poly_C1_V1[round][i], poly_C1_V1[round][i], tmp);
            }
        }
        num_terms >>= 1;
        update_V(C_1, num_terms, r[log_num + logN + 1 - round - 1]);
        update_V(tau_N, num_terms, r[log_num + logN + 1 - round - 1]);
        update_V(tau_2N, num_terms, r[log_num + logN + 1 - round - 1]);
    }
    mpz_set(C1r, C_1[0]);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("4.3. Proof generated.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);


    // Vf verifies the proof. (4-4)
    st = clock();
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r[i + 1 + logN + 1]);
    stat = stat ? sum_check_verification(rir, poly_C1_V1, V1z, 4, log_num + logN + 1, tmp_r, "GKR V1->C1") : 0;
    if (stat) {
        printf("....GKR V1->C1 sumcheck verified.\n");
        for (int i = 0; i < log_num; i ++)
            mpz_set(tmp_r[i], r[i + logN + 1 + 1]);
        evaluate_beta(br, z, tmp_r, log_num);
        evaluate_tau2(tr, val, r, logN + 1);
        mod_mult(tmp, C1r, br);
        mod_mult(tmp, tmp, tr);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V1->C1 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        }
    }

    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r[i + 1 + logN + 1]);
    mod_add(V2z, poly_C1_V2[0][0], poly_C1_V2[0][1]);
    stat = stat ? sum_check_verification(rir, poly_C1_V2, V2z, 4, log_num + logN + 1, tmp_r, "GKR V2->C1") : 0;
    if (stat) {
        printf("....GKR V2->C1 sumcheck verified.\n");
        evaluate_tau(t2r, val, r, logN + 1);
        mod_mult(tmp, C1r, br);
        mod_mult(tmp, tmp, t2r);
        if (mpz_cmp(rir, tmp)) {
            printf("GKR V2->C1 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
                mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        }
    }

    for (int i = 0; i < log_num; i ++) 
        mpz_set(tmp_r[i], r[i + 1 + logN + 1]);
    stat = stat ? sum_check_verification(rir, poly_V3_V2, V2z, 4, log_num, tmp_r, "GKR V2->V3") : 0;
    if (stat) {
        printf("....GKR V2->V3 sumcheck verified.\n");
        mod_mult(tmp, br, V_3[0]);
        mod_mult(tmp, tmp, V_3[1]);
        if (mpz_cmp(tmp, rir)) {
            printf("GKR V2->V3 Fail: One point reduction: %s %s.\n",
                mpz_get_str(0, digit_rep, rir),
               mpz_get_str(0, digit_rep, tmp));
            stat = 0;
        }
    }

    if (stat) {
        mod_1neg(rir, r[logN + 1]);
        mod_mult(rir, rir, V_3[0]);
        mod_mult(tmp, r[logN + 1], V_3[1]);
        mod_add(rir, rir, tmp); // rir = (1-r0)V3(R0)+r0V3(R1) (= V3(Rr))
        printf("..GKR all verified.\n..Returns C1r (%s) and rir (%s).\n",
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, rir));
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.4. Proof verified.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);



    // Vf checks the consistencies: rir and V3r / C1r and commit (5)

    // Vf checks whether C1r is consistent with the commit. (5-1)
    st = clock();
    kxi_eval(Cr_evalpts, CMnum, &(tmp_Cr[log_ck]));
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;

    // Open commit
    commit_open(CK_out, CK_in, Cr_evalpts, commits, CKnum, CMnum, sks);

    st = clock();
	evaluate_V(COMMIT_R, CK_out, log_ck, tmp_rCr);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;

    ////////// Commit
    st = clock();
    if(mpz_cmp(COMMIT_R, C1r)) {
        printf("Commit Fail: Inconsistent with the committed value: %s %s\n", 
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, COMMIT_R));
        stat = 0;
    } else {
        printf("....Consistent with the committed value.\n");
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.1. Check consistency with the commit.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);
    
    // Vf checks whether rir is consistent with the precomputed V3r. (5-2).
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    } else {
        printf("....Consistent with the original input.\n");
    }

    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5.2. Check consistency with the input layer.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    if (stat) {
        printf("..Consistency verified.\n");
    }
   
    
    printf("Result: %d\n", stat);
    printf("Total time:\nVf: %f ms.\nPv: %f ms.\n", TIME_VERIFIER, TIME_PROVER);
    return stat;
}


int main_FX_mult(int argc, char **argv)
{
    clock_t st;
    double elapsed_time;
    if(argc != 4) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), num = 1 << log_num, stat = 1;
    init_field(bits_of_prime, logN_in); //init field

    mpz_t val, V2z, V3r, tmp, rir;
    mpz_inits(val, V2z, V3r, tmp, rir, 0);
    mpz_t *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp3 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(tmp1[i], tmp2[i], tmp3[i], 0);

    mpz_t *beta = (mpz_t *) malloc(sizeof(mpz_t) * num);          
    for (int i = 0; i < num; i ++)
        mpz_init(beta[i]);
    mpz_t **poly_V3_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * log_num);
    for (int i = 0; i < log_num; i ++) {
        poly_V3_V2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly_V3_V2[i][j], 0);
    }
    mpz_t b[4], V30[4], V31[4];
    for (int i = 0; i < 4; i ++)
        mpz_inits(b[i], V30[i], V31[i], 0);

    mpz_t **v_3 = (mpz_t **) malloc(sizeof(mpz_t *) * num * 2),
          **v_2 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          *V_3 = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *V_2 = (mpz_t *) malloc(sizeof(mpz_t) * num);
    for (int i = 0 ; i < num; i ++) {
        mpz_inits(V_2[i], 0);
        v_2[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N * 2; j ++)
            mpz_init(v_2[i][j]);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_3[i]);
        v_3[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_3[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_3 be the array of polynomials of degree N. (1)
    }
    printf("1.1. Memory allocated and Input generated.\n");

    // Pv evaluate the circuit over polynomial rings. (2)
    // Assume that Vf gives the **coefficients** of input polynomials; hence, Pv should do fft.
    st = clock();
    for (int i = 0; i < num; i ++) {
        fft(tmp1, 2 * N, v_3[2 * i], N);
        fft(tmp2, 2 * N, v_3[2 * i + 1], N);
        fourier_mult(tmp3, tmp1, tmp2, 2 * N);
        ifft(v_2[i], tmp3, 2 * N);
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("2. Circuit evaluated.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3-1)
    // Vf also should do the same procedure but only on the input and output layer. (3-2)
    st = clock();
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N);
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    TIME_VERIFIER += elapsed_time;
    printf("3. Constructed OPR circuit.\n  Both +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Vf fix the randomness, and precompute the desired values of input and output layer. (4-1)
    mpz_t *z = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *tmp_r = (mpz_t *) malloc(sizeof(mpz_t) * log_num);
    for (int i = 0; i < log_num; i ++)
        mod_init_random(z[i]);
    for (int i = 0; i < log_num + 1; i ++)
        mod_init_random(r[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_init_set(tmp_r[i], r[i + 1]);
    printf("4.1. Randomness fixed.\n"); 

    // Vf evaluate the MLE at the chosenm points. (4-2)
    st = clock();
    evaluate_V(V2z, V_2, log_num, z);
    evaluate_V(V3r, V_3, log_num + 1, r);
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.2. Verifier precomputes the correct answer.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // The GKR Protocol. (4-3)
    // Pv generates the proof for the circuit. 
    // Pv gives poly_V3_V2[][], V_3[0], V_3[1] to Vf.
    st = clock();
    initialize_beta(beta, log_num, z);
    uint64 num_terms = 1ULL << log_num;
    for (int round = 0; round < log_num; round ++) {
        for (uint64 p = 0; p < num_terms / 2; p ++) {
            mlmap_evaluation_N(b, 4, num_terms, p, beta);
            mlmap_evaluation_N(V30, 4, num_terms << 1, p << 1, V_3);
            mlmap_evaluation_N(V31, 4, num_terms << 1, (p << 1) + 1, V_3);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b[i], V30[i]);
                mod_mult(tmp, tmp, V31[i]);
                mod_add(poly_V3_V2[round][i], poly_V3_V2[round][i], tmp);
            }
        }
        num_terms >>= 1;
        update_V(V_3, num_terms << 1, r[log_num + 1- round - 1]);
        update_V(beta, num_terms, r[log_num + 1- round - 1]);
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_PROVER += elapsed_time;
    printf("4.3. Proof generated.\n  Prover +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Vf verifies the proof. (4-4)
    st = clock();
    stat = stat ? sum_check_verification(rir, poly_V3_V2, V2z, 4, log_num, tmp_r, "GKR") : 0;
    if (stat) {
        printf("....GKR sumcheck verified.\n");
        evaluate_beta(tmp, z, tmp_r, log_num);
        mod_mult(tmp, tmp, V_3[0]);
        mod_mult(tmp, tmp, V_3[1]); // tmp = beta(z;r)V3(r0)V3(r1)
    }

    if (mpz_cmp(tmp, rir)) {
        printf("GKR Fail: One point reduction: %s %s.\n",
            mpz_get_str(0, digit_rep, rir),
            mpz_get_str(0, digit_rep, tmp));
    	stat = 0;
    }

    if (stat) {
        printf("....GKR one point reduction verified.\n");
        mod_1neg(rir, r[0]);
        mod_mult(rir, rir, V_3[0]);
        mod_mult(tmp, r[0], V_3[1]);
        mod_add(rir, rir, tmp); // rir = (1-r0)V3(R0)+r0V3(R1) (= V3(Rr))
    }

    if (stat)
        printf("..GKR All Verified\n");
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("4.4. Proof verified.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    // Check the consistency between the result of GKR and the original input. (5)
    st = clock();
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    }
    if (stat) {
        printf("..Consistent with the original input.\n");
    }
    elapsed_time = (double) ((double) clock() - st) / (CLOCKS_PER_SEC) * 1000;
    TIME_VERIFIER += elapsed_time;
    printf("5. Check consistency with the input.\n  Verifier +%f\n  %f/%f.\n-----------------------------------\n", elapsed_time, TIME_VERIFIER, TIME_PROVER);

    printf("Result: %d\n", stat);
    printf("Total time:\nVf: %f ms.\nPv: %f ms.\n", TIME_VERIFIER, TIME_PROVER);
    return stat;
}

