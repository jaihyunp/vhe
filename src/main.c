#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include "mlmap.h"
#include "vheaan.h"
#include <stdio.h>
#include <stdlib.h>

int main_rescale(int argc, char **argv);
int main_mult(int argc, char **argv);

int main(int argc, char **argv)
{
    return main_rescale(argc, argv);
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
          **v_0d = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          **v_0r = (mpz_t **) malloc(sizeof(mpz_t *) * (num * 2)),
          *V_1c = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2 * N * ubits * 4)),
          *V_0c = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2 * 2 * N * ubits * 4));
    for (int i = 0; i < num; i ++) {
        v_1r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N); //In order to mult, this have to store 2N points
        for (int k = 0; k < (int) (2 * N); k ++) 
            mpz_init(v_1r[i][k]);
        v_1d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
        for (int k = 0; k < (int) (2 * N); k ++)
            mpz_init(v_1d[i][k]);
        for (int j = 0; j < 4; j ++) {
            v_2[i * 4 + j] = (mpz_t *) malloc(sizeof(mpz_t) * N);
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
        }
    }


    //mle construction
    mpz_t P, *evk1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N), *evk2 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N);
    for (int i = 0; i < (int) N; i ++)
        mpz_inits(evk1[i], evk2[i], 0);
    mpz_init(P);

    mpz_set_ui(P, 2);
    mod_pow_ui(P, P, bits);
    mod_sub_ui(mask, P, 1);

    mod_pow(MAX, P, P);
    for (int i = 0; i < (int) N; i ++) {
        mpz_urandomm(tmp1[i], STATE, MAX);
        mpz_urandomm(tmp2[i], STATE, MAX); //Compute 2N points of N deg poly
    }
    fft(evk1, 2 * N, tmp1, N);
    fft(evk2, 2 * N, tmp2, N); //evk <- N deg poly w/ coeff of 2B bits

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
                    mpz_set_ui(V_1c[(i * 2 * N + j) * ubits * 4 + k], 1);
                } else {
                    mpz_set_ui(V_1c[(i * 2 * N + j) * ubits * 4 + k], 0);
                }
            }
        } //set V_1c
        for (int j = 0; j < (int) N; j ++) {           
            mod_sub(tmp2[j], tmp1[j], tmp1[j + N]);
            mpz_and(tmp2[j], tmp2[j], mask);
        }
        fft(v_1r[i], 2 * N, tmp2, N); //set v_1r to be reduced v_1d

        for (int j = 0; j < (int) N; j ++)
            mod_mult(tmp1[j], v_1l[2 * i][j], P);
        fourier_mult(tmp2, v_1r[i], evk1, 2 * N);
        fourier_add(v_0d[2 * i], tmp1, tmp2, 2 * N);

        for (int j = 0; j < (int) N; j ++)
            mod_mult(tmp1[j], v_1l[2 * i + 1][j], P);
        fourier_mult(tmp2, v_1r[i], evk2, 2 * N);
        fourier_add(v_0d[2 * i + 1], tmp1, tmp2, 2 * N); //v_0d <- P*v1l + v_1r*evk

        for (int ii = 0; ii < 2; ii ++) {
            ifft(tmp1, v_0d[2 * i + ii], 2 * N);
            for (int j = 0; j < (int) N; j ++) {
                for (int k = 0; k < ubits * 4 ; k ++) {
                    if (mpz_tstbit(tmp1[j], k)) {
                        mpz_set_ui(V_0c[((i * 2 + ii) * N + j)* ubits * 4 + k], 1);
                    } else {
                        mpz_set_ui(V_0c[((i * 2 + ii) * N + j)* ubits * 4 + k], 0);
                    }
                }
            } // set V_0c
            for (int j = 0; j < (int) N; j ++) {
                mod_sub(tmp2[j], tmp1[j], tmp1[j + N]);
                mpz_tdiv_q(tmp2[j], tmp2[j], P);
                mpz_and(tmp2[j], tmp2[j], mask);
            }
            fft(v_0r[2 * i + ii], 2 * N, tmp2, N); //set v_0r to be reduced v_1d
        }
   }


    //mle construction
    mpz_t *V_2 = (mpz_t *) malloc(sizeof(mpz_t) * (num * 4)),
          *V_1l = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          *V_1r = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_1d = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_0d = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          *V_0r = (mpz_t *) malloc(sizeof(mpz_t) * (num * 2)),
          EVK1, EVK2;
    mpz_inits(EVK1, EVK2, 0);
    for (int i = 0; i < num; i ++) {
        mpz_inits(V_1d[i], V_1r[i], 0);
        for (int j = 0; j < 4; j ++)
            mpz_init(V_2[i * 4 + j]);
        for (int j = 0; j < 2; j ++)
            mpz_inits(V_1l[2 * i + j], V_0d[2 * i + j], V_0r[2 * i + j]);
    }
    for (int i = 0; i < num; i ++) {
        for (int j = 0; j < 4; j ++)
            fourier_extrapolate(V_2[4 * i + j], val, v_2[4 * i + j], 2 * N);
        for (int j = 0; j < 4; j ++) {
            fourier_extrapolate(V_1l[2 * i + j], val, v_1l[2 * i + j], 2 * N);
            fourier_extrapolate(V_0d[2 * i + j], val, v_0d[2 * i + j], 2 * N);
            fourier_extrapolate(V_0r[2 * i + j], val, v_0r[2 * i + j], 2 * N);
        }
        fourier_extrapolate(V_1r[i], val, v_1r[i], 2 * N);
        fourier_extrapolate(V_1d[i], val, v_1d[i], 2 * N);
    }
    fourier_extrapolate(EVK1, val, evk1, 2 * N);
    fourier_extrapolate(EVK2, val, evk2, 2 * N);
   
    
    //randomness
    mpz_t *r_0b = (mpz_t *) malloc(sizeof(mpz_t) * ()),
          *r_0c = (mpz_t *) malloc(sizeof(mpz_t) * ()),
          *r_1b = (mpz_t *) malloc(sizeof(mpz_t) * ()),
          *r_1c = (mpz_t *) malloc(sizeof(mpz_t) * ()),
          *r_2 = (mpz_t *) malloc(sizeof(mpz_t) * ());

    mpz_t *r_2 = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2)),
          *r_1l = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *r_1r = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *r_1d = (mpz_t *) malloc(sizeof(mpz_t) * log_num),
          *r_0d = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *r_0r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1)),
          *r_1c = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + log_bits)),
          *r_0c = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2 + logN + log_bits)),
          *r_1b = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + log_bits)),
          *r_0b = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2 + logN + log_bits));
    for (int i = 0; i < log_num; i ++) {
        mod_init_random(r_1r[i]);
        mod_init_random(r_1d[i]);
    }
    for (int i = 0; i < log_num + 1; i ++) {
        mod_init_random(r_1l[i]);
        mod_init_random(r_0d[i]);
        mod_init_random(r_0r[i]);
    }
    for (int i = 0; i < log_num + 2; i ++)
        mod_init_random(r_2[i]);
    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++)
        mod_init_random(r_1c[i]);
    for (int i = 0; i < (int) (log_num + 2 + logN + log_bits); i ++)
        mod_init_random(r_0c[i]);

    //gkr
    stat = gkr_cipher_mult(
        V_0r, V_0c, V_0d,       V_1l, V_1r, V_1c, V_1d,       V_2, 
        r_0r, r_0c, r_0d, r_0b, r_1l, r_1r, r_1c, r_1d, r_1b, r_2, 
        P, EVK1, EVK2, 
        log_num, bits, val
    );


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
