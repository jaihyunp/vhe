#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include "mlmap.h"
#include "vheaan.h"
#include <stdio.h>
#include <stdlib.h>

int main_rescale(int argc, char **argv);
int main_mult(int argc, char **argv);
int main_poly_mult(int argc, char **argv);
int main(int argc, char **argv)
{
    return main_poly_mult(argc, argv);
}

int main_poly_mult(int argc, char **argv)
{
    if(argc != 4) {
        printf("Sample usage: ./a.out [bit of prime (8k-1)] [logN] [log_num].\n");
        return 0;
    }
    long bits_of_prime = atol(argv[1]);
    int logN_in = atoi(argv[2]), log_num = atoi(argv[3]), num = 1 << log_num, stat = 1;
    init_field(bits_of_prime, logN_in); //init field

    mpz_t val, V1z, V2r, tmp;
    mpz_inits(val, V1z, V2r, tmp, 0);
    mpz_t *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp3 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(tmp1[i], tmp2[i], tmp3[i], 0);


    mpz_t **v_2 = (mpz_t **) malloc(sizeof(mpz_t *) * num * 2),
          **v_1 = (mpz_t **) malloc(sizeof(mpz_t *) * num),
          *V_2 = (mpz_t *) malloc(sizeof(mpz_t *) * num * 2),
          *V_2_0 = (mpz_t *) malloc(sizeof(mpz_t *) * num),
          *V_2_1 = (mpz_t *) malloc(sizeof(mpz_t *) * num),
          *V_1 = (mpz_t *) malloc(sizeof(mpz_t *) * num);
    for (int i = 0 ; i < num; i ++) {
        mpz_inits(V_1[i], V_2_0[i], V_2_1[i], 0);
        v_1[i] = (mpz_t *) malloc(sizeof(mpz_t) * N * 2);
        for (uint64 j = 0; j < N * 2; j ++)
            mpz_init(v_1[i][j]);
    }
    for (int i = 0 ; i < num * 2; i ++) {
        mpz_init(V_2[i]);
        v_2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
        for (uint64 j = 0; j < N; j ++)
            mod_init_random(v_2[i][j]);
        for (uint64 j = N; j < N * 2; j ++)
            mpz_init_set_ui(v_2[i][j], 0);  // v_2 be the array of polynomials of degree N.
    }
    printf("Memory allocated\n");

    // Assume that Vf gives the **coefficients** of input polynomials; hence, Pv should execute fft.
    // Note that Pv does NOT need to ifft, since he can do extrapolate on random points with the fft values.
    // This would be done by Pv.
    for (int i = 0; i < num; i ++) {
        fft(tmp1, 2 * N, v_2[2 * i], N);
        fft(tmp2, 2 * N, v_2[2 * i + 1], N);
        fourier_mult(tmp3, tmp1, tmp2, 2 * N);
        ifft(v_1[i], tmp3, 2 * N);
        printf(" %d aa\n", i);
    } 
printf("aa\n");
    // Vf takes a random point val from our finite field, and Pv construct a reduced MLE circuit with val.
    // THis would be done by both Pv and Vf.
    mod_random(val);
    for (int i = 0; i < num; i ++)
        coeff_evaluate(V_1[i], val, v_1[i], 2 * N);
    for (int i = 0; i < num * 2; i ++)
        coeff_evaluate(V_2[i], val, v_2[i], N);
printf("aa\n");

    // Vf fix the randomness, and precompute the desired values of input and output layer.
    // This would be done by Vf.
    mpz_t *z = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *tmp_r = (mpz_t *) malloc(sizeof(mpz_t) * log_num);
    for (int i = 0; i < log_num; i ++)
        mod_init_random(z[i]);
    for (int i = 0; i < log_num + 1; i ++)
        mod_init_random(r[i]);
    for (int i = 0; i < log_num; i ++)
        mpz_init_set(tmp_r[i], r[i + 1]);
    evaluate_V(V1z, V_1, log_num, z);
    evaluate_V(V2r, V_2, log_num + 1, r);

    // The GKR Protocol.
    mpz_t *beta = (mpz_t *) malloc(sizeof(mpz_t) * num);          
    for (int i = 0; i < num; i ++)
        mpz_init(beta[i]);
    initialize_beta(beta, log_num, z);
    mpz_t **poly = (mpz_t **) malloc(sizeof(mpz_t *) * log_num);
    for (int i = 0; i < log_num; i ++) {
        poly[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly[i][j], 0);
    }
    mpz_t b[4], V20[4], V21[4];
    for (int i = 0; i < 4; i ++)
        mpz_inits(b[i], V20[i], V21[i], 0);

    //Proof
    uint64 num_terms = 1ULL << log_num;
    for (int round = 0; round < log_num; round ++) {
        for (uint64 p = 0; p < num_terms / 2; p ++) {
            mlmap_evaluation_N(b, 4, num_terms, p, beta);
            mlmap_evaluation_N(V20, 4, num_terms << 1, p << 1, V_2);
            mlmap_evaluation_N(V21, 4, num_terms << 1, (p << 1) + 1, V_2);

            for (int i = 0; i < 4; i ++) {
                mod_mult(tmp, b[i], V20[i]);
                mod_mult(tmp, tmp, V21[i]);
                mod_add(poly[round][i], poly[round][i], tmp);
            }
        }
        num_terms >>= 1;
        update_V(V_2, num_terms << 1, r[log_num + 1- round - 1]);
        update_V(beta, num_terms, r[log_num + 1- round - 1]);
    }

    //Prover now gives proof (consists of poly, V_2[0], V_2[1])

    //Verify
    mpz_t rir;
    mpz_inits(rir, 0);
    stat = stat ? sum_check_verification(rir, poly, V1z, 4, log_num, tmp_r, "Poly Mult") : 0;
    if (stat)
        printf("Mult over Poly Verified\n");

    evaluate_beta(tmp, z, tmp_r, log_num);
    mod_mult(tmp, tmp, V_2[0]);
    mod_mult(tmp, tmp, V_2[1]);

    if(mpz_cmp(tmp, rir))
        printf("Fail GKR Protocol: the last round %s %s\n",
            mpz_get_str(0, digit_rep, rir),
            mpz_get_str(0, digit_rep, tmp));

    mod_1neg(rir, r[0]);
    mod_mult(rir, rir, V_2[0]);
    mod_mult(tmp, r[0], V_2[1]);
    mod_add(rir, rir, tmp);

    if(mpz_cmp(V2r, rir))
        printf("Fail: Input does not match: %s %s\n",
            mpz_get_str(0, digit_rep, V2r),
            mpz_get_str(0, digit_rep, rir));

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
    mpz_t *poly = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
          *tmp2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) (2 * N); i ++)
        mpz_inits(poly[i], tmp1[i], tmp2[i], 0);

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
    mpz_t *poly = (mpz_t *) malloc(sizeof(mpz_t) * N);
    for (int i = 0; i < (int) N; i ++)
        mpz_init(poly[i]);


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
            mpz_set(poly[j], v_d[i * N + j]);
        coeff_evaluate(V_d[i], val, poly, N);
        for (int j = 0; j < (int) N; j ++)
            mpz_set(poly[j], v_r[i * N + j]);
        coeff_evaluate(V_r[i], val, poly, N);
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
        mpz_clear(poly[i]);
    for (int i = 0; i < num * 2; i ++)
        mpz_clears(V_d[i], V_r[i], NULL);
    free(poly);
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
