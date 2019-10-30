#include "field.h"
#include "vheaan.h"
/* Main for mult
*/
int main(int argc, char** argv)
{

    if(argc != 2) {
        printf("./a.out [logN]\n");
        exit(0);
    }
  
    clock_t t = clock();
    int stat = 0;
    int logN = atoi(argv[1]), N = 1 << logN;
    mpz_t tmp1, tmp2, tmp3, val, ri, ri0;
    mpz_inits(tmp1, tmp2, tmp3, val, ri, ri0, NULL);

    mpz_init(PRIME);
    mpz_set_str(PRIME, "205311158772819412817284986118022354041", 10);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));

    mpz_set_ui(val, rand());
    mod(val, val);

    mpz_t z[1], r[2];
    mpz_inits(z[0], r[0], r[1], NULL);

    mpz_set_ui(z[0], rand());
    mod(z[0], z[0]);
    mpz_set(r[0], z[0]);
    mpz_set_ui(r[1], rand());
    mod(r[1], r[1]);

    //run through entire Muggles protocol with prover
    mpz_t* V1 = (mpz_t*) malloc(sizeof(mpz_t) << 2);
    mpz_t* V0 = (mpz_t*) malloc(sizeof(mpz_t) << 1);

    for(int i = 0; i < 4; i ++) {
        mpz_init_set_ui(V1[i], 0);
        if(i < 2)
            mpz_init_set_ui(V0[i], 0);
    }
        
    for(int i = 0; i < 4; i ++) {
        mpz_set_ui(V1[i], rand());
        mod(V1[i], V1[i]);
        if(i >= 2) {
            mod_add(V0[i - 2], V1[i], V1[i - 2]);
        }
        
    }

    printf("%f sec on generating dataset\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
  
    evaluate_V(ri, V0, 1, z);
    evaluate_V(ri0, V1, 2, r);
    printf("%f sec on preprocess on in/out layers\n", (double)(clock() - t)/CLOCKS_PER_SEC);

    stat = sum_check_cipher_add(ri, V1, ri, r[1], z[0]);
    
    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
        else 
            printf("Success\n");
    } else {
        printf("Failed inside the layer\n");
    }

    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
  
    return 1;
}

/* Main for cmult
int main(int argc, char** argv)
{

    if(argc != 2) {
        printf("./a.out [logN]\n");
        exit(0);
    }
  
    clock_t t = clock();
    int stat = 0;
    int logN = atoi(argv[1]), N = 1 << logN;
    mpz_t tmp1, tmp2, tmp3, val, ri, ri0;
    mpz_inits(tmp1, tmp2, tmp3, val, ri, ri0, NULL);

    mpz_init(PRIME);
    mpz_set_str(PRIME, "205311158772819412817284986118022354041", 10);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));

    mpz_set_ui(val, rand());
    mod(val, val);

    mpz_t z[1];
    mpz_inits(z[0], NULL);

    mpz_set_ui(z[0], rand());
    mod(z[0], z[0]);
    

    //run through entire Muggles protocol with prover
    mpz_t* V1 = (mpz_t*) malloc(sizeof(mpz_t) << 1);
    mpz_t* V0 = (mpz_t*) malloc(sizeof(mpz_t) << 1);

    for(int i = 0; i < 2; i ++) {
        mpz_init_set_ui(V0[i], 0);
        mpz_init_set_ui(V1[i], 0);
    }
        
    for(int i = 0; i < 2; i ++) {
        mpz_set_ui(V1[i], rand());
        mod(V1[i], V1[i]);
        mod_mult(V0[i], V1[i], val);
//        printf("V1[%d] %s\nV0[%d] %s\n", i, mpz_get_str(NULL, 10, V1[i]), i, mpz_get_str(NULL, 10, V0[i]));
    }

    printf("%f sec on generating dataset\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
  
    evaluate_V(ri, V0, 1, z);
    evaluate_V(ri0, V1, 1, z);
    printf("%f sec on preprocess on in/out layers\n", (double)(clock() - t)/CLOCKS_PER_SEC);

    stat = sum_check_cipher_cmult(ri, V1, val, ri, z[0]);
//int sum_check_cipher_cmult(mpz_t rop, mpz_t *V, const mpz_t c, const mpz_t ri, const mpz_t z);
    
    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
        else 
            printf("Success\n");
    } else {
        printf("Failed inside the layer\n");
    }

    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
  
    return 1;
}
*/
/* Main for adding
int main(int argc, char** argv)
{

    if(argc != 2) {
        printf("./a.out [logN]\n");
        exit(0);
    }
  
    clock_t t = clock();
    int stat = 0;
    int logN = atoi(argv[1]), N = 1 << logN;
    mpz_t tmp1, tmp2, tmp3, val, ri, ri0;
    mpz_inits(tmp1, tmp2, tmp3, val, ri, ri0, NULL);

    mpz_init(PRIME);
    mpz_set_str(PRIME, "205311158772819412817284986118022354041", 10);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));

    mpz_set_ui(val, rand());
    mod(val, val);

    mpz_t z[1], r[2];
    mpz_inits(z[0], r[0], r[1], NULL);

    mpz_set_ui(z[0], rand());
    mod(z[0], z[0]);
    mpz_set(r[0], z[0]);
    mpz_set_ui(r[1], rand());
    mod(r[1], r[1]);

    //run through entire Muggles protocol with prover
    mpz_t* V1 = (mpz_t*) malloc(sizeof(mpz_t) << 2);
    mpz_t* V0 = (mpz_t*) malloc(sizeof(mpz_t) << 1);

    for(int i = 0; i < 4; i ++) {
        mpz_init_set_ui(V1[i], 0);
        if(i < 2)
            mpz_init_set_ui(V0[i], 0);
    }
        
    for(int i = 0; i < 4; i ++) {
        mpz_set_ui(V1[i], rand());
        mod(V1[i], V1[i]);
        if(i >= 2) {
            mod_add(V0[i - 2], V1[i], V1[i - 2]);
        }
        
    }

    printf("%f sec on generating dataset\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
  
    evaluate_V(ri, V0, 1, z);
    evaluate_V(ri0, V1, 2, r);
    printf("%f sec on preprocess on in/out layers\n", (double)(clock() - t)/CLOCKS_PER_SEC);

    stat = sum_check_cipher_add(ri, V1, ri, r[1], z[0]);
    
    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
        else 
            printf("Success\n");
    } else {
        printf("Failed inside the layer\n");
    }

    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
  
    return 1;
}
*/

/* Main for rescaling
int main(int argc, char **argv)
{
    if((argc != 4)) {
        printf("./a.out [log degree] [bits before] [bits after]\n");
        return 0;
    } 
 
    clock_t t = clock();
    int stat = 0;
    int logN = atoi(argv[1]), N = 1 << logN, d1 = atoi(argv[2]), d2 = atoi(argv[3]), lb = 0;
    while((d1 - 1) >> ++ lb);
    int b = 1 << lb, l = 1 + logN + lb;
    mpz_t tmp1, tmp2, tmp3, val, ri, ri0;
    mpz_inits(tmp1, tmp2, tmp3, val, ri, ri0, NULL);

    mpz_init(PRIME);
    mpz_set_str(PRIME, "205311158772819412817284986118022354041", 10);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));

    mpz_set_ui(val, rand());
    mod(val, val);


    mpz_t* rn = (mpz_t*) malloc(1 * sizeof(mpz_t));
    mpz_t* rd = (mpz_t*) malloc(logN * sizeof(mpz_t));
    mpz_t* rb = (mpz_t*) malloc(lb * sizeof(mpz_t));

    mpz_t* zn = (mpz_t*) malloc(1 * sizeof(mpz_t));
    mpz_t* zdb = (mpz_t*) malloc((logN + lb) * sizeof(mpz_t));

    for(int i = 0; i < 1; i ++) {
        mpz_init_set_ui(rn[i], rand());
        mod(rn[i], rn[i]);
    }


    for(int i = 0; i < logN; i ++) {
        mpz_init_set_ui(rd[i], rand());
        mod(rd[i], rd[i]);
    }


    for(int i = 0; i < lb; i ++) {
        mpz_init_set_ui(rb[i], rand());
        mod(rb[i], rb[i]);
    }

    for(int i = 0; i < 1; i ++) {
        mpz_init_set_ui(zn[i], rand());
        mod(zn[i], zn[i]);
    }
    
     for(int i = 0; i < logN + lb; i ++) {
        mpz_init_set_ui(zdb[i], rand());
        mod(zdb[i], zdb[i]);
    }


    mpz_t* V1 = (mpz_t*) malloc(sizeof(mpz_t) << 1);
    mpz_t* V0 = (mpz_t*) malloc(sizeof(mpz_t) << 1);
    mpz_t* coeffs1 = (mpz_t*) malloc(sizeof(mpz_t) << logN);
    mpz_t* coeffs2 = (mpz_t*) malloc(sizeof(mpz_t) << logN);
    mpz_t* C = (mpz_t*) malloc(sizeof(mpz_t) << (1 + logN + lb));    

    for(int i = 0; i < 2; i ++) {
        mpz_init_set_ui(V1[i], 0);
        mpz_init_set_ui(V0[i], 0);
    }
    for(int i = 0; i < N; i ++) {
        mpz_inits(coeffs1[i], coeffs2[i], NULL);
    }
    for(int i = 0; i < 2 * N * b; i ++) {
        mpz_init_set_ui(C[i], 0);
    }
        

    for(int i = 0; i < 2; i ++) {
        mpz_set_ui(tmp2, 1);
        for(int j = 0; j < N; j ++) {
            mpz_set_ui(coeffs1[j], 0);
            mpz_set_ui(coeffs2[j], 0);
            mpz_set_ui(tmp1, 1);
            mpz_set_ui(tmp3, 1);
            for(int k = 0; k < d1; k ++) {
                if(rand() & 1) {
                    mpz_set_ui(C[i * N * b + j * b + k], 1);
                    mpz_add(coeffs1[j], coeffs1[j], tmp1);

                    if(k >= (d1 - d2)) {
                        mpz_add(coeffs2[j], coeffs2[j], tmp3);
                    }
                }
                mod_add(tmp1, tmp1, tmp1);
                if(k >= (d1 - d2))
                    mod_add(tmp3, tmp3, tmp3);
            }

            mod_mult(tmp1, coeffs1[j], tmp2);
            mod_add(V1[i], V1[i], tmp1);
            
            mod_mult(tmp1, coeffs2[j], tmp2);
            mod_add(V0[i], V0[i], tmp1);
            
            mod_mult(tmp2, tmp2, val);
        }
    }

    printf("%f sec on generating dataset\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
  
    evaluate_V(ri, V0, 1, zn);
    evaluate_V(ri0, V1, 1, zn);
    printf("%f sec on preprocess on in/out layers\n", (double)(clock() - t)/CLOCKS_PER_SEC);

    stat = sum_check_cipher_rescale(ri, C, ri, logN, d1, d2, rn, rd, rb, val, zn, zdb);

    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
    } else {
        printf("Failed inside the layer\n");
    }
    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
} 
*/
