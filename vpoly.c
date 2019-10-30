/*************************************
//Jai Hyun Park
//October 27, 2019.
//Implementation of polynomial rounding circuit using commitments
**************************************/

#include "vheaan.h"
int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);


    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check %s %s %s\n", str, mpz_get_str(NULL, 10, tmp), mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, PRIME));
        stat = 0;
    }
    
    if(stat) {
        extrapolate(extrap_val, poly[0], degree, r[num_rounds - 1]);
        for(int round = 1; round < num_rounds; round ++) {
            mod_add(tmp, poly[round][0], poly[round][1]);
            if(mpz_cmp(tmp, extrap_val)) {
                printf("Fail sumcheck %s :: %dth round %s %s\n", str, round, mpz_get_str(NULL, 10, tmp), mpz_get_str(NULL, 10, extrap_val));
                stat = 0;
                break;
            } else {
                extrapolate(extrap_val, poly[round], degree, r[num_rounds - 1 - round]);
            }
        }
    }

    if(stat) {
        mpz_set(rop, extrap_val);
    } else {
        mpz_set_ui(rop, ERROR);
    }

    mpz_clears(tmp, extrap_val, NULL);

    return stat;
}

// log # poly = ln
// log degree = ld
// d1 bits -> d2 bits
int sum_check_poly_rounding(mpz_t rop, mpz_t *C, const mpz_t ri, const int ln, const int ld, const int d1, const int d2, const mpz_t *r_n, const mpz_t *r_d, const mpz_t *r_b, const mpz_t val, const mpz_t *z_n, const mpz_t *z_db)
{
    clock_t timer = clock();
    int stat = 1;
    int lb = 0;
    while((d1 - 1) >> ++ lb);
    int l = ln + ld + lb;
    mpz_t tmp1, tmp2, cr, a1r, a2r, b1r, b2r, tr;
    mpz_inits(tmp1, tmp2, cr, a1r, a2r, b1r, b2r, tr, NULL);
    mpz_t *r = (mpz_t*) malloc(sizeof(mpz_t) * (ln + ld + lb)),
          *z = (mpz_t*) malloc(sizeof(mpz_t) * (ln + ld + lb));
    for(int i = 0; i < lb; i ++)
        mpz_init_set(r[i], r_b[i]);
    for(int i = 0; i < ld; i ++)
        mpz_init_set(r[i + lb], r_d[i]);
    for(int i = 0; i < ln; i ++)
        mpz_init_set(r[i + lb + ld], r_n[i]);
    for(int i = 0; i < ld + lb; i ++)
        mpz_init_set(z[i], z_db[i]);
    for(int i = 0; i < ln; i ++)
        mpz_init_set(z[i + ld + lb], z_n[i]);
        
    //commit
    evaluate_V(cr, C, l, r);
    evaluate_alpha(a1r, d1, r_b, lb);
    evaluate_alpha(a2r, d1 - d2, r_b, lb);
    evaluate_beta(b1r, z_n, r_n, ln);
    evaluate_beta(b2r, z, r, l);
    evaluate_tau(tr, val, r_d, ld); 
    printf("%f sec on pre process of Verifier\n", (double)(clock() - timer) / CLOCKS_PER_SEC);
    timer = clock();
 
    
    //sum check
    
    //init
    mpz_t b1[3], b2[4], c[4], a1[3], a2[3], t[3];
    for(int i = 0; i < 3; i ++)
        mpz_inits(b1[i], b2[i], c[i], a1[i], a2[i], t[i], NULL);
    mpz_inits(b2[3], c[3], NULL);
    
    mpz_t *beta1 = (mpz_t*) malloc(sizeof(mpz_t) << ln),
          *beta2 = (mpz_t*) malloc(sizeof(mpz_t) << l),
          *tau = (mpz_t*) malloc(sizeof(mpz_t) << ld),
          *alpha1 = (mpz_t*) malloc(sizeof(mpz_t) << lb),
          *alpha2 = (mpz_t*) malloc(sizeof(mpz_t) << lb);
    for(uint64 i = 0; i < (1ULL << ln); i ++)
        mpz_init(beta1[i]);
    for(uint64 i = 0; i < (1ULL << l); i ++)
        mpz_init(beta2[i]);
    for(uint64 i = 0; i < (1ULL << ld); i ++)
        mpz_init(tau[i]);
    for(uint64 i = 0; i < (1ULL << lb); i ++)
        mpz_inits(alpha1[i], alpha2[i], NULL);

    initialize_beta(beta1, ln, z_n);
    initialize_beta(beta2, l, z);
    initialize_tau(tau, ld, val);
    initialize_alpha(alpha1, lb, d1);
    initialize_alpha(alpha2, lb, d1 - d2);

    mpz_t **poly_b = (mpz_t**) malloc(sizeof(mpz_t*) * l),
          **poly_d = (mpz_t**) malloc(sizeof(mpz_t*) * l),
          **poly_r = (mpz_t**) malloc(sizeof(mpz_t*) * l);
    for(int i = 0; i < l; i ++){
        poly_b[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        poly_d[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        poly_r[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        for(int j = 0; j < 4; j ++){
            mpz_init_set_ui(poly_b[i][j], 0);
        }
        for(int j = 0; j < 3; j ++){
            mpz_init_set_ui(poly_d[i][j], 0);
            mpz_init_set_ui(poly_r[i][j], 0);
        }
    }
    printf("%f sec on init of Prover\n", (double)(clock() - timer) / CLOCKS_PER_SEC);
    timer = clock();

    //proof
    uint64 num_terms = 1ULL << l;
    for(int round = 0; round < l; round ++) {
        
        for(uint64 j = 0; j < num_terms / 2; j ++) {
            
            if(round < ln) {
                linear_evaluation(b1, 3, num_terms >> (lb + ld), j >> (lb + ld), beta1);
                linear_evaluation(t, 3, 0, (j >> lb) & ((1 << ld) - 1), tau);
                linear_evaluation(a1, 3, 0, j & ((1 << lb) - 1), alpha1);
                linear_evaluation(a2, 3, 0, j & ((1 << lb) - 1), alpha2);
            } else if (round < ln + ld) {
                linear_evaluation(b1, 3, 0, 0, beta1);
                linear_evaluation(t, 3, num_terms >> lb, j >> lb, tau);
                linear_evaluation(a1, 3, 0, j & ((1 << lb) - 1), alpha1);
                linear_evaluation(a2, 3, 0, j & ((1 << lb) - 1), alpha2); 
            } else { 
                linear_evaluation(b1, 3, 0, 0, beta1);
                linear_evaluation(t, 3, 0, 0, tau);
                linear_evaluation(a1, 3, num_terms, j, alpha1);
                linear_evaluation(a2, 3, num_terms, j, alpha2);
            }

            linear_evaluation(b2, 4, num_terms, j, beta2);
            linear_evaluation(c, 4, num_terms, j, C);
            
            for(int k = 0; k < 3; k ++) {
                mod_mult(tmp1, b1[k], c[k]);
                mod_mult(tmp1, t[k], tmp1);
                mod_mult(tmp2, tmp1, a2[k]);//tmp2<-bk*tk*ck*a2k
                mod_mult(tmp1, tmp1, a1[k]);//tmp1<-bk*tk*ck*a1k

                mod_add(poly_d[round][k], poly_d[round][k], tmp1);
                mod_add(poly_r[round][k], poly_r[round][k], tmp2);
            }

            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b2[k]); //tmp1<-bk*ck*(ck-1)
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }
        }
        
        num_terms >>= 1;

        if(round < ln) {
            update_V(beta1, num_terms >> (lb + ld), r[l - 1 - round]);
        } else if (round < ln + ld) {
            update_V(tau, num_terms >> lb, r[l - 1 - round]);
        } else {
            update_V(alpha1, num_terms, r[l - 1 - round]);
            update_V(alpha2, num_terms, r[l - 1 - round]);
        }

        update_V(C, num_terms, r[l - 1 - round]);
        update_V(beta2, num_terms, r[l - 1 - round]);
        
    }

    printf("%f sec on generating proof of Prover\n", (double)(clock() - timer) / CLOCKS_PER_SEC);
    timer = clock();

    //verify
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, l, r, "V_b") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, b2r);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_b output does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    mod_add(tmp1, poly_d[0][0], poly_d[0][1]); 
    stat = stat ? sum_check_verification(tmp2, poly_d, tmp1, 3, l, r, "V_d") : 0;
    mod_mult(tmp1, a1r, cr);
    mod_mult(tmp1, tmp1, b1r);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_d output does not match\n");
    }

    mod_add(tmp1, poly_d[0][0], poly_d[0][1]); 
    mpz_set_ui(tmp2, 2);
    mod_pow(tmp2, tmp2, d1 - d2);
    mod_mult(tmp2, tmp2, ri);
    mod_sub(tmp1, tmp1, tmp2);
    stat = stat ? sum_check_verification(tmp2, poly_r, tmp1, 3, l, r, "V_r") : 0;
    mod_mult(tmp1, a2r, cr);
    mod_mult(tmp1, tmp1, b1r);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_r output does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    printf("%f sec on verification of Verifier\n", (double)(clock() - timer) / CLOCKS_PER_SEC);
    timer = clock();

    //return
    if(stat) {
        printf("Layer Verified\n");
        mod_add(rop, poly_d[0][0], poly_d[0][1]);
    } else
        mpz_set_ui(rop, 0);

    //free
    mpz_clears(tmp1, tmp2, cr, a1r, a2r, b1r, b2r, tr, NULL);
    for(int i = 0; i < l; i ++)
        mpz_clear(r[i]);
    for(int i = 0; i < l; i ++)
        mpz_clear(z[i]);
    free(r);
    free(z);
    for(int i = 0; i < 3; i ++)
        mpz_clears(b1[i], b2[i], c[i], a1[i], a2[i], t[i], NULL);
    mpz_clears(b2[3], c[3], NULL);
    
    for(uint64 i = 0; i < (1ULL << ln); i ++)
        mpz_clear(beta1[i]);
    for(uint64 i = 0; i < (1ULL << l); i ++)
        mpz_clear(beta2[i]);
    for(uint64 i = 0; i < (1ULL << ld); i ++)
        mpz_clear(tau[i]);
    for(uint64 i = 0; i < (1ULL << lb); i ++)
        mpz_clears(alpha1[i], alpha2[i], NULL);
    free(beta1);
    free(beta2);
    free(alpha1);
    free(alpha2);
    free(tau);

    for(int i = 0; i < l; i ++) {
        for(int j = 0; j < 3; j ++)
            mpz_clears(poly_b[i][j], poly_d[i][j], poly_r[i][j], NULL);
        mpz_clear(poly_b[i][3]);
        free(poly_b[i]);
        free(poly_d[i]);
        free(poly_r[i]);
    }
    free(poly_b);
    free(poly_d);
    free(poly_r);

    return stat;
}


int sum_check_rounding_fast(mpz_t rop, mpz_t *C, const mpz_t ri, const int d, const int d1, const int d2, const mpz_t *r, const mpz_t* z)
{
    clock_t t = clock(), tall = clock();
    int stat = 1;
    int l = 0;
    while((d1 - 1) >> ++ l);
    mpz_t tmp1, tmp2, cr, a1r, a2r, br, b[4], c[4], a1[3], a2[3], ri_b; 
    mpz_inits(tmp1, tmp2, cr, a1r, a2r, br, ri_b, NULL);
    for(int i = 0; i < 4; i ++)
        mpz_inits(b[i], c[i], NULL);
    for(int i = 0; i < 3; i ++)
        mpz_inits(a1[i], a2[i], NULL);

    //commit
    evaluate_V(cr, C, d + l, r);
    
    //alpha
    mpz_t *r_head = (mpz_t*) malloc(l * sizeof(mpz_t)),
         *r_tail = (mpz_t*) malloc(d * sizeof(mpz_t));
    for(int i = 0; i < l; i ++)
        mpz_init_set(r_head[i], r[i]);
    for(int i = 0; i < d; i ++)
        mpz_init_set(r_tail[i], r[i + l]);
    
    evaluate_alpha(a1r, d1, r_head, l);
    evaluate_alpha(a2r, d1 - d2, r_head, l);
    mod_inv_pow2(tmp1, d1 - d2);
    mod_sub(a2r, a1r, a2r);
    mod_mult(a2r, a2r, tmp1);    

    //beta
    evaluate_beta(br, z, r, l + d);

    printf("%f sec til prover's work\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();

    mpz_t *alpha1 = (mpz_t*) malloc(sizeof(mpz_t) << l),
          *alpha2 = (mpz_t*) malloc(sizeof(mpz_t) << l),
          *betavals = (mpz_t*) malloc(sizeof(mpz_t) << (l + d));

    for(int i = 0; i < (1 << l); i ++) {
        if(i < d1) {
            mpz_init_set_ui(alpha1[i], 2);
            mod_pow(alpha1[i], alpha1[i], i);
        } else {
            mpz_init_set_ui(alpha1[i], 0);
        }
    }
    for(int i = 0; i < (1 << l); i ++) {
        if((i >= (d1 - d2)) && (i < d1)) {
            mpz_init_set_ui(alpha2[i], 2);
            mod_pow(alpha2[i], alpha2[i], i - d1 + d2);
        } else {
            mpz_init_set_ui(alpha2[i], 0);
        }
    }
    
    for(int i = 0; i < (1 << l + d); i ++)
        mpz_init(betavals[i]);
    initialize_beta(betavals, l + d, z); 
    
    mpz_t **poly_b = (mpz_t**) malloc(sizeof(mpz_t*) * (d + l)),
          **poly_d = (mpz_t**) malloc(sizeof(mpz_t*) * l),
          **poly_r = (mpz_t**) malloc(sizeof(mpz_t*) * l);
    for(int i = 0; i < d + l; i ++) {
        poly_b[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        for(int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly_b[i][j], 0);
    }
    for(int i = 0; i < l; i ++) {
        poly_d[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        poly_r[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        for(int j = 0; j < 3; j ++){
            mpz_init_set_ui(poly_d[i][j], 0);
            mpz_init_set_ui(poly_r[i][j], 0);
        }
    }

    uint64 num_terms = 1ULL << (d + l);
    for(int round = 0; round < d; round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++) {

            linear_evaluation(b, 4, num_terms, j, betavals);
            linear_evaluation(c, 4, num_terms, j, C);

            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]); //tmp1<-bk*ck*(ck-1)
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }
        }
        num_terms >>= 1;
        update_V(betavals, num_terms, r_tail[d - 1 - round]);
        update_V(C, num_terms, r_tail[d - 1 - round]);
    }

    printf("%f sec til pro1\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
 

    //verifier do:
    mpz_set_ui(tmp1, 0);
    if(!sum_check_verification(ri_b, poly_b, tmp1, 4, d, r_tail, "inter V_b")) {
        stat = 0;
        printf("Fail :: inter V_b\n");
    }

    printf("%f sec til ver1\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
 

    //prover do:
    for(int i = 0; i < l; i ++)
        for(int j = 0; j < 4; j ++)
            mpz_set_ui(poly_b[i][j], 0);

    if(stat) {
        for(int round = 0; round < l; round ++) {
            for(uint64 j = 0; j < num_terms / 2; j ++) {

                linear_evaluation(b, 4, num_terms, j, betavals);
                linear_evaluation(c, 4, num_terms, j, C);
                linear_evaluation(a1, 3, num_terms, j, alpha1);
                linear_evaluation(a2, 3, num_terms, j, alpha2);
                
                for(int k = 0; k < 3; k ++) {
                    mod_mult(tmp1, a1[k], c[k]);
                    mod_add(poly_d[round][k], poly_d[round][k], tmp1);

                    mod_mult(tmp1, a2[k], c[k]);
                    mod_add(poly_r[round][k], poly_r[round][k], tmp1);
                }

                for(int k = 0; k < 4; k ++) {
                    mpz_sub_ui(tmp1, c[k], 1);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_mult(tmp1, tmp1, b[k]);
                    mod_add(poly_b[round][k], poly_b[round][k], tmp1);
                }
            
            }
            
            num_terms >>= 1;
            update_V(betavals, num_terms, r_head[l - 1 - round]);
            update_V(C, num_terms, r_head[l - 1 - round]);
            update_V(alpha1, num_terms, r_head[l - 1 - round]);
            update_V(alpha2, num_terms, r_head[l - 1 - round]);
        }
    }

    printf("%f sec til pro2\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
 

    //verifier do:
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, br);
    stat = stat ? sum_check_verification(tmp2, poly_b, ri_b, 4, l, r_head, "V_b") : 0;
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_b output does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    mod_mult(tmp1, a2r, cr);
    stat = stat ? sum_check_verification(tmp2, poly_r, ri, 3, l, r_head, "V_r") : 0;
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_r output does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    mod_add(tmp1, poly_d[0][0], poly_d[0][1]); 
    stat = stat ? sum_check_verification(tmp2, poly_d, tmp1, 3, l, r_head, "V_d") : 0;
    mod_mult(tmp1, a1r, cr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: V_d output does not match\n");
    }

    if(stat)
        mod_add(rop, poly_d[0][0], poly_d[0][1]);
    else
        mpz_set_ui(rop, 0);

    printf("%f sec til ver2\n", (double)(clock() - t)/CLOCKS_PER_SEC);
 

    for(int i = 0; i < 4; i ++)
        mpz_clears(b[i], c[i], NULL);
    for(int i = 0; i < 3; i ++)
        mpz_clears(a1[i], a2[i], NULL);
    for(int i = 0; i < d + l; i ++) { 
        for(int j = 0; j < 4; j ++)
            mpz_clear(poly_b[i][j]);
        free(poly_b[i]);
    }
    free(poly_b);
    for(int i = 0; i < l; i ++) {
        for(int j = 0; j < 3; j ++)
            mpz_clears(poly_d[i][j], poly_r[i][j], NULL);
        free(poly_d[i]);
        free(poly_r[i]);
    }
    free(poly_d);
    free(poly_r);

    for(int i = 0; i < (1 << l + d); i ++)
        mpz_clear(betavals[i]);
    for(int i = 0; i < (1 << l); i ++)
        mpz_clears(alpha1[i], alpha2[i], NULL);
    free(alpha1);
    free(alpha2);
    free(betavals);
    for(int i = 0; i < l; i ++)
        mpz_clear(r_head[i]);
    for(int i = 0; i < d; i ++)
        mpz_clear(r_tail[i]);
    free(r_head);
    free(r_tail);
    mpz_clears(tmp1, tmp2, cr, a1r, a2r, br, ri_b, NULL);

    printf("%f sec for a layer\n", (double)(clock() - tall)/CLOCKS_PER_SEC);
    return stat;
}
