#include "typedefs.h"
#include "field.h"
#include "parameters.h"
#include "polynomial.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    if(argc != 4) {
        printf("./a.out [Prime bit] [logN_in] [mode]\n");
        printf("mode 0: ROU ~\n");
        printf("mode 1: FFT ~ \n");
        printf("mode 2: ADD ~ \n");
        printf("mode 3: MULT ~\n");
        printf("mode 4: EXTRAP ~\n");
        return 0;
    }
    
    /* Init */
    int logN_in = atoi(argv[2]);
    long prime_bit = atol(argv[1]);
    int mode = atoi(argv[3]);
    init_field(prime_bit, logN_in);

    mpz_t *v1, *v2, *v, *c1, *c2, *c;
    v = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    v1 = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    v2 = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    c = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    c1 = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    c2 = (mpz_t*) malloc(sizeof(mpz_t) * N * 2);
    for(uint64 i = 0; i < 2 * N; i ++)
        mpz_inits(v1[i], v2[i], v[i], c1[i], c2[i], c[i], NULL);

    mpz_t tmp, val1, val2;
    mpz_inits(tmp, val1, val2, NULL);


    if(mode <= 0) {
    
        /* Test for ROU */
        printf("====================================================\n");
        printf("Test for ROU: 4N th root of unity\n");
        printf("----------------------------------------------------\n");
        mod_pow_ui(tmp, ROU[1], 2 * N);
        if(!mpz_cmp_ui(tmp, 1))
            printf("Fail to ROU1\n");
        mod_mult(tmp, tmp, tmp);
        if(mpz_cmp_ui(tmp, 1))
            printf("Fail to ROU2\n");
        printf("ROU Test done\n");   
        printf("====================================================\n\n\n");
    }

    if(mode <= 1) {

        /* Test for FFT & iFFT */
        printf("====================================================\n");
        printf("Test for FFT & iFFT\n");
        printf("----------------------------------------------------\n");
        for(uint64 i = 0; i < N; i ++) {
            mod_random(c[i]);
        }
        fft(v, N, c, N);
        printf(" fft done\n");
        for(uint64 i = 0; i < N; i ++) {
            coeff_evaluate(tmp, ROU[i * 4], c, N);
            if(mpz_cmp(tmp, v[i]))
                printf("Fail to FFT\n");
        }
        printf("FFT Test done\n");
        
        ifft(c1, v, N);
        printf(" ifft done\n");
        for(uint64 i = 0; i < N; i ++) {
            if(mpz_cmp(c1[i], c[i]))
                printf("Fail to iFFT\n");
        }
        printf("iFFT Test done\n");
        printf("====================================================\n\n\n");
    }

    if(mode <= 2) {

        /* Test for Add */
        printf("====================================================\n");
        printf("Test for Add\n");
        printf("----------------------------------------------------\n");
        for(uint64 i = 0; i < N; i ++) {
            mod_random(c1[i]);
            mod_random(c2[i]);
        }

        fft(v1, N, c1, N);
        fft(v2, N, c2, N);
        printf(" fft done\n");
        fourier_add(v, v1, v2, N);    
        printf("f-add done\n");
        coeff_add(c, c1, c2, N);
        printf("c-add done\n");

        if(polynomial_cmp(c, N, v, N))
            printf("Fail to add\n");
        printf("Poly Add Test done\n");
        printf("====================================================\n\n\n");
    }

    if(mode <= 3) {

        /* Test for Mult */
        printf("====================================================\n");
        printf("Test for Mult\n");
        printf("----------------------------------------------------\n");
        for(uint64 i = 0; i < N; i ++) {
            mod_random(c1[i]);
            mod_random(c2[i]);
        }

        fft(v1, 2 * N, c1, N);
        fft(v2, 2 * N, c2, N);
        printf(" fft done\n");
        fourier_mult(v, v1, v2, 2 * N);
        printf("f-mult done\n");
        coeff_mult(c, c1, c2, N);
        printf("c-mult done\n");
        
        if(polynomial_cmp(c, 2 * N, v, 2 * N))
            printf("Fail to mult\n");
        printf("Poly Mult Test done\n");
        printf("====================================================\n\n\n");
    }

    if(mode <= 4) {

        /* Test for Eval & Extrap */
        printf("====================================================\n");
        printf("Test for Extrap\n");
        printf("----------------------------------------------------\n");
        for(uint64 i = 0; i < 2 * N; i ++) {
            mod_random(c1[i]);
            mod_random(c2[i]);
        }

        fft(v1, N, c1, N);
        fft(v2, 2 * N, c2, 2 * N);
        printf(" fft done\n");

        mod_random(tmp);
        
        coeff_evaluate(val1, tmp, c1, N);
        printf(" c1 done\n");
        fourier_extrapolate(val2, tmp, v1, N);
        printf(" f1 done\n");
        if(mpz_cmp(val1, val2))
            printf("Fail to extrap N deg\n");
        printf("N degree done\n");

        coeff_evaluate(val1, tmp, c2, 2 * N);
        printf(" c2 done\n");
        fourier_extrapolate(val2, tmp, v2, 2 * N);
        printf(" f2 done\n");
        if(mpz_cmp(val1, val2))
            printf("Fail to extrap 2N deg\n");
        printf("2N degree done\n");

        printf("Fourier Extrap Test done\n");
        printf("====================================================\n\n\n");
    }

    /* Free */
    mpz_clears(tmp, val1, val2, NULL);
    for(uint64 i = 0; i < N * 2; i ++){
        mpz_clears(v1[i], v2[i], v[i], c1[i], c2[i], c[i], NULL);
    }
    free(v1);
    free(v2);
    free(v);
    free(c1);
    free(c2);
    free(c);
    free_field();
}

