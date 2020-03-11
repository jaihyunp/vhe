#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include "mlmap.h"
#include "vc.h"
#include <stdio.h>
#include <stdlib.h>

int main_FX_mult(int argc, char **argv);
int main_FXoverPhi_mult(int argc, char **argv);
int main_RXoverPhi_mult(int argc, char **argv);

int main(int argc, char **argv)
{
//    return main_FX_mult(argc, argv);
//    return main_FXoverPhi_mult(argc, argv);
    return main_RXoverPhi_mult(argc, argv);
}


int main_RXoverPhi_mult(int argc, char **argv)
{
    if(argc != 5) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num] [bits].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), num = 1 << log_num, stat = 1, log_bits = 0, bits = atoi(argv[4]), ubits;
    while ((bits - 1) >> ++ log_bits);
    ubits = 1 << log_bits;

    init_field(bits_of_prime, logN_in); //init field

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
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_3[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_3 be the array of polynomials of degree N. (1)
    }
    for (int i = 0; i < num * (int) N * 2; i ++)
        mpz_init(C_1[i]);
    for (int i = 0; i < num * (int) N * ubits * 4; i ++)
        mpz_init(C_0[i]);
    mod_ui_pow_ui(BITMAX, 2, 4 * ubits - 1);
    mod_ui_pow_ui(Q, 2, bits);
    printf("..Memory allocated.\n");


    // Pv evaluates the circuit and commit the required values. (2)
    // Pv evaluates the circuit over polynomial rings. (2-1)
    // Assume that Vf gives and takes the **coefficients** of input polynomials; hence, Pv should do (i)fft.
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
    //Commit (2-2)

    //////////
    // TODO //
    ////////// Commit

    printf("..Circuit evaluated.\n"); 


    
    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3)
    // Vf also should do the same procedure but only on the input and output layer.
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_0[i], val, v_0[i], N); // Both Pv and Vf
    for (int i = 0; i < num; i ++)
        for (uint64 j = 0; j < N; j ++)
            mpz_set(V_1[i * N + j], v_1[i][j]);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N); // Pv
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N); // Both Pv and Vf
    printf("..MLE constructed.\n"); 


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
    printf("....Randomness fixed.\n"); 

    //////////
    // TODO //
    ////////// Commit
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r1[i + 1]);
    evaluate_V(COMMIT_R_1, C_1, log_num + logN + 1, tmp_r);

    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(tmp_r[i], r3[i]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(tmp_r[i + log_bits + 2], r2[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + log_bits + 2], r1[i + 1]);
    evaluate_V(COMMIT_R_0, C_0, log_num + logN + log_bits + 2, tmp_r);


    // Vf evaluates the MLE of input and output layers at the chosen point. (4-2)
    evaluate_V(V0z, V_0, log_num, z1);
    evaluate_V(V3r, V_3, log_num + 1, r1);
    printf("....Prover Precomputed.\n"); 


    // Pv generates the proof for the circuit. (4-3) 
    // Pv gives poly_V3_V2[][][], V_3[0], V_3[1], and C1r to Vf.
    initialize_beta(beta1, log_num, z1);
    initialize_beta(beta2, logN, z2);
    initialize_beta(beta3, log_bits + 2, z3);
    initialize_alpha(alpha_4Q, log_bits + 2, bits * 4);
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
    for (int i = 0; i < 4; i ++)
        mod_set_ui(sign[i], 1 - 2 * i);
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
    printf("....GKR Proof generated.\n");


    // Vf verifies the proof. (4-4)
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
            printf("GKR V1->bits Fail: One point reduction: %s %s.\n",
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
        evaluate_alpha(a4qr, 4 * bits, r3, log_bits + 2);
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
    mod_add(V1z, poly_C1_V1[0][0], poly_C1_V1[0][1]);
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
        printf("..GKR all verified. Returns C0r (%s), C1r (%s) and rir (%s).\n",
                mpz_get_str(0, digit_rep, C0r),
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, rir));
    }



    // Vf checks the consistencies: rir and V3r / C1r, C0r and commit (5)
    // Vf checks whether C*r are consistent with the commits. (5-1)

    //////////
    // TODO //
    //////////
    if(mpz_cmp(COMMIT_R_0, C0r)) {
        printf("Commit Fail: Inconsistent with the committed value: %s %s\n", 
                mpz_get_str(0, digit_rep, C0r),
                mpz_get_str(0, digit_rep, COMMIT_R_0));
        stat = 0;
    } else {
        printf("....Consistent with the committed value.\n");
    }
    if(mpz_cmp(COMMIT_R_1, C1r)) {
        printf("Commit Fail: Inconsistent with the committed value: %s %s\n", 
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, COMMIT_R_1));
        stat = 0;
    } else {
        printf("....Consistent with the committed value.\n");
    }

    
    // Vf checks whether rir is consistent with the precomputed V3r. (5-2).
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    } else {
        printf("....Consistent with the original input.\n");
    }

    if (stat) {
        printf("..Consistency verified.\n");
    }
 
    printf("Result: %d\n", stat);
    return stat;
}


//Implementation of Verifiable Multiplication over Z_p[X]/(X^N + 1), where p is prime and N is power of two.
int main_FXoverPhi_mult(int argc, char **argv)
{
    if(argc != 4) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), num = 1 << log_num, stat = 1;
    init_field(bits_of_prime, logN_in); //init field

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

    printf("..Memory allocated.\n");


    // Pv evaluates the circuit and commit the required values. (2)
    // Pv evaluates the circuit over polynomial rings. (2-1)
    // Assume that Vf gives and takes the **coefficients** of input polynomials; hence, Pv should do (i)fft.
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
    //Commit (2-2)

    //////////
    // TODO //
    //////////

    printf("..Circuit evaluated.\n"); 


    
    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3)
    // Vf also should do the same procedure but only on the input and output layer.
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_1[i], val, v_1[i], N); // Both Pv and Vf
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N); // Pv
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N); // Both Pv and Vf
    printf("..MLE constructed.\n"); 


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
    printf("....Randomness fixed.\n"); 

    //////////
    // TODO //
    //////////
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(tmp_r[i], r[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_set(tmp_r[i + logN + 1], r[i + 1 + logN + 1]);
    evaluate_V(COMMIT_R, C_1, log_num + logN + 1, tmp_r);

    // Vf evaluates the MLE of input and output layers at the chosen point. (4-2)
    evaluate_V(V1z, V_1, log_num, z);
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(tmp_r[i], r[i + logN + 1]);
    evaluate_V(V3r, V_3, log_num + 1, tmp_r);
    printf("....Prover Precomputed.\n"); 

    // Pv generates the proof for the circuit. (4-3) 
    // Pv gives poly_V3_V2[][][], V_3[0], V_3[1], and C1r to Vf.
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
    printf("....GKR Proof generated.\n");

    // Vf verifies the proof. (4-4)
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
        printf("..GKR all verified. Returns C1r (%s) and rir (%s).\n",
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, rir));
    }



    // Vf checks the consistencies: rir and V3r / C1r and commit (5)
    // Vf checks whether C1r is consistent with the commit. (5-1)

    //////////
    // TODO //
    //////////
    if(mpz_cmp(COMMIT_R, C1r)) {
        printf("Commit Fail: Inconsistent with the committed value: %s %s\n", 
                mpz_get_str(0, digit_rep, C1r),
                mpz_get_str(0, digit_rep, COMMIT_R));
        stat = 0;
    } else {
        printf("....Consistent with the committed value.\n");
    }
    
    // Vf checks whether rir is consistent with the precomputed V3r. (5-2).
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    } else {
        printf("....Consistent with the original input.\n");
    }

    if (stat) {
        printf("..Consistency verified.\n");
    }
    
    
    
    printf("Result: %d\n", stat);
    return stat;
}


int main_FX_mult(int argc, char **argv)
{
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
    printf("..Memory allocated.\n");

    // Pv evaluate the circuit over polynomial rings. (2)
    // Assume that Vf gives the **coefficients** of input polynomials; hence, Pv should do fft.
    for (int i = 0; i < num; i ++) {
        fft(tmp1, 2 * N, v_3[2 * i], N);
        fft(tmp2, 2 * N, v_3[2 * i + 1], N);
        fourier_mult(tmp3, tmp1, tmp2, 2 * N);
        ifft(v_2[i], tmp3, 2 * N);
    }
    printf("..Circuit evaluated.\n"); 

    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val. (3-1)
    // Vf also should do the same procedure but only on the input and output layer. (3-2)
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], 2 * N);
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_3[i], val, v_3[i], N);
    printf("..MLE constructed.\n"); 

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
    printf("....Randomness fixed.\n"); 

    // Vf evaluate the MLE at the chosenm points. (4-2)
    evaluate_V(V2z, V_2, log_num, z);
    evaluate_V(V3r, V_3, log_num + 1, r);
    printf("....Prover Precomputed.\n"); 

    // The GKR Protocol. (4-3)
    // Pv generates the proof for the circuit. 
    // Pv gives poly_V3_V2[][], V_3[0], V_3[1] to Vf.
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
    printf("....GKR Proof generated.\n");

    // Vf verifies the proof. (4-4)
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

    // Check the consistency between the result of GKR and the original input. (5)
    if (mpz_cmp(V3r, rir)) {
        printf("GKR Fail: Inconsistent with the input: %s %s.\n",
            mpz_get_str(0, digit_rep, V3r),
            mpz_get_str(0, digit_rep, rir));
    	stat = 0;
    }
    if (stat) {
        printf("..Consistent with the original input.\n");
    }

    printf("Result: %d\n", stat);
    return stat;
}

