#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include "mlmap.h"
#include "vheaan.h"
#include <stdio.h>
#include <stdlib.h>

int main_rescale(int argc, char **argv);
int main_mult(int argc, char **argv);
int main_FX_mult(int argc, char **argv);
int main_FXoverPhi_mult(int argc, char **argv);
int main_RXoverPhi_mult(int argc, char **argv);
int main(int argc, char **argv)
{
    return main_FXoverPhi_mult(argc, argv);
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

    mpz_t val, V1z, V2z, V3r, C1r, tmp, rir, t0N, br, tr, t2r, COMMIT_R_1, COMMIT_R_0;
    mpz_inits(val, V1z, V2z, V3r, C1r, tmp, rir, t0N, br, tr, t2r, COMMIT_R_1, COMMIT_R_0, 0);
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
          **poly_C1_V2 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1)),
          **poly_C0_V0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_bits + 2)),
          **poly_C0_V1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2));
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
            mpz_init(v_2[i][j]);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_3[i]);
        v_3[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_3[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_2 be the array of polynomials of degree N. (1)
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
    evaluate_V(COMMIT_R_1, C_1, log_num + logN + 1, tmp_r);

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
    mod_pow_ui(t0N, val, N);
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
    if(mpz_cmp(COMMIT_R_0, C1r)) {
        printf("Commit Fail: Inconsistent with the committed value: %s %s\n", 
                mpz_get_str(0, digit_rep, COMMIT_R_0),
                mpz_get_str(0, digit_rep, COMMIT_R_0));
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

    mpz_t val, V1z, V2z, V3r, C1r, tmp, rir, t0N, br, tr, t2r, COMMIT_R;
    mpz_inits(val, V1z, V2z, V3r, C1r, tmp, rir, t0N, br, tr, t2r, COMMIT_R, 0);
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
            mpz_init(v_2[i][j]);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_3[i]);
        v_3[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_3[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_3[i][j], 0);  // v_2 be the array of polynomials of degree N. (1)
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
    mod_pow_ui(t0N, val, N);
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
            mpz_init_set_ui(v_3[i][j], 0);  // v_2 be the array of polynomials of degree N. (1)
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

//for mult;
int main_mult(int argc, char **argv)
{
    if (argc != 5) {
        printf("Usage: ./a.out [bit of prime (8k-1)] [logN] [log_num] [bits].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), bits = atoi(argv[4]), log_bits = 0, num = 1 << log_num, stat = 1;
    init_field(bits_of_prime, logN_in); //init field
    while ((bits - 1) >> ++ log_bits);
    int ubits = 1 << log_bits;
    mpz_t base1, base2, val, MAX, mask;
    mpz_inits(base1, base2, val, MAX, mask, NULL);
    mpz_t *poly_V3_V2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(poly_V3_V2[i], tmp1[i], tmp2[i], 0);

    //set val
    mpz_urandomm(val, STATE, PRIME);

    //circ evaluation & construction
    mpz_t **v_2 = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 4)),
          **v_1l = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          **v_1r = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_1d = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_1s = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          **v_0r = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          **v_0d = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          **v_0s = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          *V_1c0 = (mpz_t *) malloc(sizeof(mpz_t) * (num * N * 2 * ubits * 4)),
          *V_1c1 = (mpz_t *) malloc(sizeof(mpz_t) * (num * N * ubits * 4)),
          *V_0c0 = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2 * N * 2 * ubits * 4)),
          *V_0c1 = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2 * N * ubits * 4));

    for (int i = 0; i < (int) (num * 2 * N * ubits * 4); i ++) 
        mpz_init(V_1c0[i]);
    for (int i = 0; i < (int) (num * N * ubits * 4); i ++) 
        mpz_init(V_1c1[i]);
    for (int i = 0; i < (int) (num * 2 * N * 2 * ubits * 4); i ++) 
        mpz_init(V_0c0[i]);
    for (int i = 0; i < (int) (num * 2 * N * ubits * 4); i ++) 
        mpz_init(V_0c1[i]);

    for (int i = 0; i < num; i ++) { //In order to mult, we need 2N points for all poly
        v_1r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
        for (int k = 0; k < (int) (2 * N); k ++) 
            mpz_init(v_1r[i][k]);
        v_1d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
        for (int k = 0; k < (int) (2 * N); k ++)
            mpz_init(v_1d[i][k]);
        v_1s[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
        for (int k = 0; k < (int) (2 * N); k ++)
            mpz_init(v_1s[i][k]);
        for (int j = 0; j < 4; j ++) {
            v_2[i * 4 + j] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
            for (int k = 0; k < (int) (2 * N); k ++)
                mpz_init(v_2[i * 4 + j][k]);
        }
        for (int j = 0; j < 2; j ++) {
            v_1l[i * 2 + j] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
            for (int k = 0; k < (int) (2 * N); k ++) 
                mpz_init(v_1l[i * 2 + j][k]);
            v_0d[i * 2 + j] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
            for (int k = 0; k < (int) (2 * N); k ++) 
                mpz_init(v_0d[i * 2 + j][k]);
            v_0r[i * 2 + j] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
            for (int k = 0; k < (int) (2 * N); k ++) 
                mpz_init(v_0r[i * 2 + j][k]);
            v_0s[i * 2 + j] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
            for (int k = 0; k < (int) (2 * N); k ++) 
                mpz_init(v_0s[i * 2 + j][k]);
        }
    }
    printf("GKR Mult: Memory Allocated\n");


    //mle construction
    mpz_t P, *evk0 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N), *evk1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(evk0[i], evk1[i], 0);
    mpz_init(P);

    mpz_set_ui(P, 2);
    mod_pow_ui(P, P, bits); //P = 2^B
    mod_sub_ui(mask, P, 1); //mask = 2^B - 1

    mod_mult(MAX, P, P); //MAX = 2^(2B)
    for (int i = 0; i < (int) N; i ++) {
        mpz_urandomm(tmp1[i], STATE, MAX);
        mpz_urandomm(tmp2[i], STATE, MAX); //Compute 2N points of N deg poly
    }
    fft(evk0, 2 * N, tmp1, N);
    fft(evk1, 2 * N, tmp2, N); //evk <- N deg poly w/ coeff of 2B bits


//    mod_mult(MAX, MAX, MAX); //MAX = 2^(4B)
    mod_ui_pow_ui(MAX, 2, 4 * ubits - 1);
    printf("P: %s\n", mpz_get_str(0, digit_rep, P));
    printf("MASK: %s\n", mpz_get_str(0, digit_rep, mask));
    printf("MAX: %s\n", mpz_get_str(0, digit_rep, MAX));
    for (int i = 0; i < num; i ++) {
        for (int j = 0; j < 4; j ++) {
            for (int k = 0; k < (int) N; k ++)
                mpz_urandomm(tmp1[k], STATE, P);
            fft(v_2[4 * i + j], 2 * N, tmp1, N); //Compute 2N points of N deg poly
        } //v_2 <- N deg poly w/ coeff of B bits = (a1,b1,a2,b2)

        fourier_mult(v_1l[2 * i], v_2[4 * i], v_2[4 * i + 2], 2 * N);
        fourier_mult(tmp1, v_2[4 * i + 1], v_2[4 * i + 2], 2 * N);
        fourier_mult(tmp2, v_2[4 * i], v_2[4 * i + 3], 2 * N);
        fourier_add(v_1l[2 * i + 1], tmp1, tmp2, 2 * N); //v_1l = (a1*a2,a1*b2+a2*b1)

        fourier_mult(v_1d[i], v_2[4 * i + 1], v_2[4 * i + 3], 2 * N); //v_1d = b1*b2
        ifft(tmp1, v_1d[i], 2 * N);
        for (int j = 0; j < (int) (2 * N); j ++) {
            for (int k = 0; k < ubits * 4; k ++) {
                if (mpz_tstbit(tmp1[j], k)) {
                    mpz_set_ui(V_1c0[(i * 2 * N + j) * ubits * 4 + k], 1);
                } else {
                    mpz_set_ui(V_1c0[(i * 2 * N + j) * ubits * 4 + k], 0);
                }
            }
        } //set V_1c0
        for (int j = 0; j < (int) N; j ++) {
            mod_sub(tmp2[j], tmp1[j], tmp1[j + N]);
            mod_add(tmp2[j], tmp2[j], MAX); //substract mod 2^(4B-1)
        }
        fft(v_1s[i], 2 * N, tmp2, N); //set v_1s
        for (int j = 0; j < (int) N; j ++) {
            for (int k = 0; k < ubits * 4; k ++) {
                if (mpz_tstbit(tmp2[j], k)) {
                    mpz_set_ui(V_1c1[(i * N + j) * ubits * 4 + k], 1);
                } else {
                    mpz_set_ui(V_1c1[(i * N + j) * ubits * 4 + k], 0);
                }
            }
        } //set V_1c1
        for (int j = 0;j < (int) N; j ++) {
            mpz_and(tmp2[j], tmp2[j], mask); //extract the last B bits
        }
        fft(v_1r[i], 2 * N, tmp2, N); //set v_1r to be reduced v_1d

        for (int j = 0; j < (int) (2 * N); j ++)
            mod_mult(tmp1[j], v_1l[2 * i][j], P);
        fourier_mult(tmp2, v_1r[i], evk0, 2 * N);
        fourier_add(v_0d[2 * i], tmp1, tmp2, 2 * N);

        for (int j = 0; j < (int) (2 * N); j ++)
            mod_mult(tmp1[j], v_1l[2 * i + 1][j], P);
        fourier_mult(tmp2, v_1r[i], evk1, 2 * N);
        fourier_add(v_0d[2 * i + 1], tmp1, tmp2, 2 * N); //v_0d <- P*v1l + v_1r*evk

        for (int ii = 0; ii < 2; ii ++) {
            ifft(tmp1, v_0d[2 * i + ii], 2 * N);
            for (int j = 0; j < (int) (2 * N); j ++) {
                for (int k = 0; k < ubits * 4 ; k ++) {
                    if (mpz_tstbit(tmp1[j], k)) {
                        mpz_set_ui(V_0c0[((i * 2 + ii) * 2 * N + j) * ubits * 4 + k], 1);
                    } else {
                        mpz_set_ui(V_0c0[((i * 2 + ii) * 2 * N + j) * ubits * 4 + k], 0);
                    }
                }
            } // set V_0c0
            for (int j = 0; j < (int) N; j ++) {
                mod_sub(tmp2[j], tmp1[j], tmp1[j + N]);
                mod_add(tmp2[j], tmp2[j], MAX); //mod 2^(4B-1)
            }
            fft(v_0s[2 * i + ii], 2 * N, tmp2, N); //set v_0s
            for (int j = 0; j < (int) N; j ++) {
                for (int k = 0; k < ubits * 4 ; k ++) {
                    if (mpz_tstbit(tmp2[j], k)) {
                        mpz_set_ui(V_0c1[((i * 2 + ii) * N + j) * ubits * 4 + k], 1);
                    } else {
                        mpz_set_ui(V_0c1[((i * 2 + ii) * N + j) * ubits * 4 + k], 0);
                    }
                }
            } // set V_0c1
            for (int j = 0; j < (int) N; j ++) {
                mpz_tdiv_q(tmp2[j], tmp2[j], P);
                mpz_and(tmp2[j], tmp2[j], mask);
            }
            fft(v_0r[2 * i + ii], 2 * N, tmp2, N); //set v_0r to be reduced v_1d
        }
    }
    printf("GKR Mult: MLE evaluated\n");


    //mle construction
    mpz_t *V_2 = (mpz_t *) malloc(sizeof(mpz_t) * (num * 4)),
          *V_1l = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          *V_1r = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_1d = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_1s = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_0r = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          *V_0d = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          *V_0s = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          EVK1, EVK0;

    mpz_inits(EVK1, EVK0, 0);
    for (int i = 0; i < num; i ++) {
        mpz_inits(V_1d[i], V_1r[i], V_1s[i], 0);
        for (int j = 0; j < 4; j ++)
            mpz_init(V_2[i * 4 + j]);
        for (int j = 0; j < 2; j ++)
            mpz_inits(V_1l[2 * i + j], V_0d[2 * i + j], V_0r[2 * i + j], V_0s[2 * i + j], 0);
    }
    
    for (int i = 0; i < num; i ++) {
        for (int j = 0; j < 4; j ++)
            fourier_extrapolate(V_2[4 * i + j], val, v_2[4 * i + j], 2 * N);
        for (int j = 0; j < 2; j ++) {
            fourier_extrapolate(V_1l[2 * i + j], val, v_1l[2 * i + j], 2 * N);
            fourier_extrapolate(V_0d[2 * i + j], val, v_0d[2 * i + j], 2 * N);
            fourier_extrapolate(V_0r[2 * i + j], val, v_0r[2 * i + j], 2 * N);
            fourier_extrapolate(V_0s[2 * i + j], val, v_0s[2 * i + j], 2 * N);
        }
        fourier_extrapolate(V_1r[i], val, v_1r[i], 2 * N);
        fourier_extrapolate(V_1d[i], val, v_1d[i], 2 * N);
        fourier_extrapolate(V_1s[i], val, v_1s[i], 2 * N);
    }
    fourier_extrapolate(EVK0, val, evk0, 2 * N);
    fourier_extrapolate(EVK1, val, evk1, 2 * N);
    printf("GKR Mult: MLE (OPR) Constructed\n");
    
    

    //randomness
    mpz_t *z_num = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *z_N = (mpz_t *) malloc(sizeof(mpz_t) * (logN + 1)),
          *z_bits = (mpz_t *) malloc(sizeof(mpz_t) * (log_bits + 2)),
          *r_num = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2)),
          *r_N = (mpz_t *) malloc(sizeof(mpz_t) * (logN + 1)),
          *r_bits = (mpz_t *) malloc(sizeof(mpz_t) * (log_bits + 2));


    for (int i = 0; i < log_num + 1; i ++)
        mod_init_random(z_num[i]);
    for (int i = 0; i < log_num + 2; i ++)
        mod_init_random(r_num[i]);
    for (int i = 0; i < (int) logN + 1; i ++) {
        mod_init_random(r_N[i]);
        mod_init_random(z_N[i]);
    }
    for (int i = 0; i < log_bits + 2; i ++) {
        mod_init_random(r_bits[i]);
        mod_init_random(z_bits[i]);
    }

    stat = gkr_cipher_mult(
        V_0r, V_0s, V_0d, V_0c0, V_0c1, V_1l, V_1r, V_1s, V_1d, V_1c0, V_1c0, V_2,
        z_num, z_N, z_bits, 
        r_num, r_N, r_bits, 
        P, EVK0, EVK1, log_num, bits, val
    );
//    //gkr
//    stat = gkr_cipher_mult(
//        V_0r, V_0s, V_0d, V_0c0, V_0c1, V_1l, V_1r, V_1s, V_1d, V_1c0, V_1c0, V_2,
//        r_num, r_0, r_1, r_2,
//        P, EVK1, EVK2, 
//        log_num, bits, val
//    );

//    for (int i = 0; i < num; i ++) {
//        mod_mult(val, V_2[i * 4], V_2[i * 4 + 2]);
//        if (mpz_cmp(val, V_1l[i * 2])) {
//            printf("Invalid V1l[0] on %d: %s / %s\n", i,
//                mpz_get_str(0, digit_rep, val),
//                mpz_get_str(0, digit_rep, V_1l[i * 2]));
//        } else {
//            printf("Valid V1l[0] on %d: %s\n", i, mpz_get_str(0, digit_rep, val));
//        }
//        printf("%d th rounding:\n%s\n%s\n%s\n%s\n", i,
//            mpz_get_str(0, digit_rep, v_1s[i][0]),
//            mpz_get_str(0, digit_rep, v_1r[i][0]),
//            mpz_get_str(0, digit_rep, v_0s[i][0]),
//            mpz_get_str(0, digit_rep, v_0r[i][0]));
//
//    }

    if (stat)
        printf("Success\n");
    else
        printf("Failed\n");

    return 0;
}


//for rescale;
int main_rescale(int argc, char **argv)
{
    if (argc != 6) {
        printf("Usage: ./a.out [bit of prime (8k-1)] [logN] [log_num] [bef_bits] [aft_bits].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), bef_bits = atoi(argv[4]), aft_bits = atoi(argv[5]), log_bits = 0, dif_bits = bef_bits - aft_bits, num = 1 << log_num, stat = 1;
    init_field(bits_of_prime, logN_in); //init field
    while ((bef_bits - 1) >> ++ log_bits);
    int ubits = 1 << log_bits;
    mpz_t base1, base2, val;
    mpz_inits(base1, base2, val, NULL);
    mpz_t *poly_V3_V2 = (mpz_t *) malloc(sizeof(mpz_t) * N);
    for (int i = 0; i < (int) N; i ++)
        mpz_init(poly_V3_V2[i]);


    //set val
    mpz_urandomm(val, STATE, PRIME);

    //circ evaluation & construction
    mpz_t *v_d = (mpz_t *) malloc(sizeof(mpz_t) * (N * num * 2)),
          *v_r = (mpz_t *) malloc(sizeof(mpz_t) * (N * num * 2)),
          *v_c = (mpz_t *) malloc(sizeof(mpz_t) * (N * num * 2 * ubits));
    
    for (int i = 0; i < (int) N * num * 2; i ++) {
        mpz_init_set_ui(v_d[i], 0);
        mpz_init_set_ui(v_r[i], 0);

        for (int j = 0; j < ubits; j ++) {
            mpz_init(v_c[i * ubits + j]);

            if (j == 0)
                mpz_set_ui(base1, 1);
            if (j == dif_bits)
                mpz_set_ui(base2, 1);

            if (j < dif_bits) {
                mpz_urandomb(v_c[i * ubits + j], STATE, 1);
                mod_mult(val, v_c[i * ubits + j], base1);
                mod_add(v_d[i], v_d[i], val);
                mpz_mul_ui(base1, base1, 2);
            } else if (j < bef_bits) {
                mpz_urandomb(v_c[i * ubits + j], STATE, 1);
                mod_mult(val, v_c[i * ubits + j], base1);
                mod_add(v_d[i], v_d[i], val);
                mod_mult(val, v_c[i * ubits + j], base2);
                mod_add(v_r[i], v_r[i], val);
                mpz_mul_ui(base1, base1, 2);
                mpz_mul_ui(base2, base2, 2);
            } else {
                mpz_set_ui(v_c[i * ubits + j], 0);
            }
        }
    }
    printf("circuit evaled\n");


    //mle construction
    mpz_t *V_d = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *V_r = (mpz_t *) malloc(sizeof(mpz_t) * num * 2);
    for (int i = 0; i < num * 2; i ++) {
        mpz_inits(V_d[i], V_r[i], NULL);
        for (int j = 0; j < (int) N; j ++)
            mpz_set(poly_V3_V2[j], v_d[i * N + j]);
        coeff_evaluate(V_d[i], val, poly_V3_V2, N);
        for (int j = 0; j < (int) N; j ++)
            mpz_set(poly_V3_V2[j], v_r[i * N + j]);
        coeff_evaluate(V_r[i], val, poly_V3_V2, N);
    }
    printf("mle constructed\n");


    //randomness
    mpz_t *r_b = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + log_bits)),
          *r_c = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + log_bits));
    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++) {
        mod_init_random(r_b[i]);
        mod_init_random(r_c[i]);
    }


    //gkr
    stat = gkr_cipher_rescale(V_r, V_d, v_c, r_c, r_b, log_num, bef_bits, aft_bits, val);


    if (stat)
        printf("Success\n");
    else
        printf("Failed\n");


    //free
    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++)
        mpz_clears(r_b[i], r_c[i], NULL);
    free(r_b);
    free(r_c);
    for (int i = 0; i < (int) N * num * 2; i ++) {
        mpz_clears(v_d[i], v_r[i], NULL);
        for (int j = 0; j < ubits; j ++)
            mpz_clear(v_c[i * ubits + j]);
    }
    for (int i = 0; i <(int) N; i ++)
        mpz_clear(poly_V3_V2[i]);
    for (int i = 0; i < num * 2; i ++)
        mpz_clears(V_d[i], V_r[i], NULL);
    free(poly_V3_V2);
    free(V_d);
    free(V_r);
    free(v_d);
    free(v_r);
    free(v_c);
    mpz_clears(base1, base2, val, NULL);

    return 0;
}



/*
    // ! This description is invalid from now: 02 Mar 2020 !
    // Description of the structure of VHE Mult:
    //
    //  This is the exact implementation of cipher mult in HEaaN.
    //  But that in most cyclotomic ring based HE schemes (such 
    //  as B/GV) are also almost same, so this diagram can be applied 
    //  to the most of HE schemes.
    //
    // (In)
    // V_2       V_1c
    // |     \  / |   \
    // V_1l  V_1d V_1r V_1b
    // |          /
    // |         /
    // |  V_0c0 /  V_0c1
    // |  | | \/ /  \    \
    // |  | | /\/    \    \
    // |  | |/ /\     \    \
    // |  | / /  \
    // |  |/|/    \
    // V_0d V_0r V_0b
    //     (Out)
    //
    // Choice of Randomness:
    //  - All MLE but commit exploits (same) 'r_num' at the head of its randomness.
    //  - Each randomness of V_1c and V_0c is totally independent from others.
    //  - V_0d, V_0r, V_1l (and EVK) (and V_0b) shares a same randomness
*/ 
