#include "vheaan.h"

//Input: C1x(t) C1y(t) C2x(t) C2y(t); z
//Layer : C1x C1y C2x C2y evkx evky
//Layer : C1xC2x C1yC2x+C1xC2y C1yC2y evkx evky
//Layer : C1xC2x C1yC2x+C1xC2y P*C1yC2yevkx P*C1yC2yevky
//Layer : 
//Output: .. :) ..
int sum_check_cipher_mult
(
    mpz_t rop, mpz_t *V, mpz_t *C1, mpz_t *C2, 
    const mpz_t ri_in,
    const mpz_t P, const mpz_t *E, 
    const int lb;
    const mpz_t *z, const mpz_t *zdb
)
{
    //init
    int stat = 1, L = 1 + (logN + 1) + (lb + 2);
    uint64 num_terms;
    mpz_t tmp1, tmp2, ri, cr, b1r, a1r, a2r, a3r, t1r, t2r;
    mpz_inits(tmp1, ri, cr, NULL);
    mpz_t **poly_b, **poly_d, **poly_r, *z, *r;

    ////layer0
    //layer0 init
    int l = 1 + (logN + 1) + (lb + 2);
    for(int i = 0; i < lb + 2; i ++) {
        r[i] = r0b[i];
        z[i] = z0db[i];
    }
    for(int i = 0; i < logN + 1; i ++) {
        r[i + lb + 2] = r0d[i];
        z[i + lb + 2] = z0db[i + lb + 2];
    }
    for(int i = 0; i < 1; i ++) {
        r[i + logN + 1 + lb + 2] = r0n[i];
        z[i + logN + 1 + lb + 2] = r0n[i];
//        z[i + logN + 1 + lb + 2] = z0n[i];
    }
    for(int i = 0; i < l; i ++){
        for(int k = 0; k < 4; k ++) {
            mpz_set_ui(poly_b[i][j], 0);
            mpz_set_ui(poly_d[i][j], 0);
            mpz_set_ui(poly_r[i][j], 0);
        }
    }


    //commit
    evaluate_V(cr, C2, l, r);
    evaluate_beta(b1r, z, r, l);
    evaluate_alpha(a1r, 1 << (lb + 2), r0b, lb + 2);
    evaluate_alpha(a2r, 1 << (lb + 1), r0b, lb + 2);
    evaluate_alpha(a3r, 1 << lb, r0b, lb + 2);
    evaluate_tau(t1r, val, r0d, logN + 1);
    evaluate_tau2(t2r, val, r0d, logN + 1);

    //Prover
    //Prover init
    initialize_beta(beta1, l, z);
    initialize_alpha(alpha1, lb + 2, 1 << (lb + 2));
    initialize_alpha(alpha2, lb + 2, 1 << (lb + 1));
    initialize_alpha(alpha3, lb + 2, 1 << lb);
    initialize_tau(tau1, logN + 1, val);
    initialize_tau2(tau2, logN + 1, val);

    //Proof
    num_terms = 1ULL << l;
    for(int round = 0; round < l; round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++) {
            linear_evaluation(c, 4, num_terms, j, C2);
            linear_evaluation(b1, 4, num_terms, j, beta);

            if(round < 1) {
            } else if (round < 1 + logN + 1) {
                linear_evaluation(t1, 3, num_terms >> (lb + 2), j >> (lb + 2), tau1);
                linear_evaluation(t2, 3, num_terms >> (lb + 2), j >> (lb + 2), tau2);
                lienar_evaluation(a1, 3, 0, j & ((1 << (lb + 2)) - 1), alpha1);
                lienar_evaluation(a2, 3, 0, j & ((1 << (lb + 2)) - 1), alpha2);
                lienar_evaluation(a3, 3, 0, j & ((1 << (lb + 2)) - 1), alpht3);
            } else {
                linear_evaluation(t1, 3, 0, 0, tau1);
                linear_evaluation(t2, 3, 0, 0, tau2);
                linear_evaluation(a1, 3, num_terms, j, alpha1);
                linear_evaluation(a2, 3, num_terms, j, alpha2);
                linear_evaluation(a3, 3, num_terms, j, alpha3);
            }

            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }
            if(round > 0) {
                for(int k = 0; k < 3; k ++){
                    mod_mult(tmp1, a1[k], t1[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_d[round - 1][k], poly_d[round - 1][k], tmp1);

                    mod_sub(tmp1, a2[k], a3[k]);
                    mod_mult(tmp1, tmp1, t2[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_r[round - 1][k], poly_r[round - 1][k], tmp1);
                }
            }

        }

        num_terms >>= 1;
        update_V(C2, num_terms, r[l - round - 1]);
        update_V(beta, num_terms, r[l - round - 1]);
        if(round < 1) {
        } else if (round < 1 + logN + 1) {
            update_V(tau1, num_terms >> (lb + 2), r[l - 1 - round]);
            update_V(tau2, num_terms >> (lb + 2), r[l - 1 - round]);
        } else {
            update_V(alpha1, num_terms, r[l - 1 - round]);
            update_V(alpha2, num_terms, r[l - 1 - round]);
            update_V(alpha3, num_terms, r[l - 1 - round]);
        }
    }

    //Verifier
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, l - 1, r, "Layer0 Vb") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, b2r);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Layer0 Vb does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    mpz_set(tmp1, ri);
    mpz_set_ui(tmp2, 2);
    mod_pow(tmp2, tmp2, lb);
    mod_mult(tmp1, tmp2);
    stat = stat ? sum_check_verification(tmp2, poly_r, tmp1, 3, l - 1, r, "Layer0 Vr") : 0;
    mod_sub(tmp1, a1r, a2r);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, tau2r);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Layer0 Vr does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }
    
    mod_add(tmp1, poly_b[0][1], poly_b[0][0]);
    stat = stat ? sum_check_verification(ri, poly_d, tmp1, 3, l - 1, r, "Layer0 Vd") : 0;
    if(stat) {
        printf("Layer0 Verified\n");
    }
   

    ////layer1
    //layer1 init
    int l = 1 + (logN + 1) + (lb + 2);
    for(int i = 0; i < lb + 2; i ++) {
        r[i] = r0b[i];
        z[i] = z0db[i];
    }
    for(int i = 0; i < logN + 1; i ++) {
        r[i + lb + 2] = r0d[i];
        z[i + lb + 2] = z0db[i + lb + 2];
    }
    for(int i = 0; i < 1; i ++) {
        r[i + logN + 1 + lb + 2] = r0n[i];
        z[i + logN + 1 + lb + 2] = r0n[i];
//        z[i + logN + 1 + lb + 2] = z0n[i];
    }

    //commit
    evaluate_V(cr, C2, l, r);
    evaluate_beta(b1r, z, r, l);
    evaluate_alpha(a1r, 1 << (lb + 2), r0b, lb + 2);
    evaluate_alpha(a2r, 1 << (lb + 1), r0b, lb + 2);
    evaluate_alpha(a3r, 1 << lb, r0b, lb + 2);
    evaluate_tau(t1r, val, r0d, logN + 1);
    evaluate_tau2(t2r, val, r0d, logN + 1);

    //Prover
    //Prover init
    initialize_beta(beta1, l, z);
    initialize_alpha(alpha1, lb + 2, 1 << (lb + 2));
    initialize_alpha(alpha2, lb + 2, 1 << (lb + 1));
    initialize_alpha(alpha3, lb + 2, 1 << lb);
    initialize_tau(tau1, logN + 1, val);
    initialize_tau2(tau2, logN + 1, val);

    //Proof
    num_terms = 1ULL << l;
    for(int round = 0; round < l; round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++) {
            linear_evaluation(c, 4, num_terms, j, C2);
            linear_evaluation(b1, 4, num_terms, j, beta);

            if(round < 1) {
            } else if (round < 1 + logN + 1) {
                linear_evaluation(t1, 3, num_terms >> (lb + 2), j >> (lb + 2), tau1);
                linear_evaluation(t2, 3, num_terms >> (lb + 2), j >> (lb + 2), tau2);
                lienar_evaluation(a1, 3, 0, j & ((1 << (lb + 2)) - 1), alpha1);
                lienar_evaluation(a2, 3, 0, j & ((1 << (lb + 2)) - 1), alpha2);
                lienar_evaluation(a3, 3, 0, j & ((1 << (lb + 2)) - 1), alpht3);
            } else {
                linear_evaluation(t1, 3, 0, 0, tau1);
                linear_evaluation(t2, 3, 0, 0, tau2);
                linear_evaluation(a1, 3, num_terms, j, alpha1);
                linear_evaluation(a2, 3, num_terms, j, alpha2);
                linear_evaluation(a3, 3, num_terms, j, alpha3);
            }

            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }
            if(round > 0) {
                for(int k = 0; k < 3; k ++){
                    mod_mult(tmp1, a1[k], t1[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_d[round - 1][k], poly_d[round - 1][k], tmp1);

                    mod_sub(tmp1, a2[k], a3[k]);
                    mod_mult(tmp1, tmp1, t2[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_r[round - 1][k], poly_r[round - 1][k], tmp1);
                }
            }

        }

        num_terms >>= 1;
        update_V(C2, num_terms, r[l - round - 1]);
        update_V(beta, num_terms, r[l - round - 1]);
        if(round < 1) {
        } else if (round < 1 + logN + 1) {
            update_V(tau1, num_terms >> (lb + 2), r[l - 1 - round]);
            update_V(tau2, num_terms >> (lb + 2), r[l - 1 - round]);
        } else {
            update_V(alpha1, num_terms, r[l - 1 - round]);
            update_V(alpha2, num_terms, r[l - 1 - round]);
            update_V(alpha3, num_terms, r[l - 1 - round]);
        }
    }

    //Verifier
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, l - 1, r, "Layer0 Vb") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, b2r);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Layer0 Vb does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }

    mpz_set(tmp1, ri);
    mpz_set_ui(tmp2, 2);
    mod_pow(tmp2, tmp2, lb);
    mod_mult(tmp1, tmp2);
    stat = stat ? sum_check_verification(tmp2, poly_r, tmp1, 3, l - 1, r, "Layer0 Vr") : 0;
    mod_sub(tmp1, a1r, a2r);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, tau2r);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Layer0 Vr does not match %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, tmp2));
    }
    
    mod_add(tmp1, poly_b[0][1], poly_b[0][0]);
    stat = stat ? sum_check_verification(ri, poly_d, tmp1, 3, l - 1, r, "Layer0 Vd") : 0;
    if(stat) {
        printf("Layer0 Verified\n");
    }
   



    //free
    mpz_clears(tmp1, NULL); 
}


//Input: C1x(t) C1y(t); z
//Params: z, c
//Output: cC1x(t) cC1y(t); z
int sum_check_cipher_cmult(mpz_t rop, mpz_t *V, const mpz_t c, const mpz_t ri, const mpz_t z)
{
    int stat = 1;
    mpz_t inv, tmp1, pf;
    mpz_inits(inv, tmp1, pf, NULL);

    //Prover
    if(mpz_cmp_ui(c, 0)){
        mod_1neg(tmp1, z);
        mod_mult(V[0], V[0], tmp1);
        mod_mult(V[1], V[1], z);
        mod_add(pf, V[0], V[1]);
    } else {
        mod_inv(inv, c);
        mod_mult(pf, ri, c);
    }

    //Verifier
    if(!mpz_cmp_ui(c, 0)) {
        if(mpz_cmp_ui(ri, 0)) {
            stat = 0;
            printf("Fail :: zero divisor\n");
        }
    } else {
        mod_mult(tmp1, pf, c);
        if(mpz_cmp(ri, tmp1)) {
            stat = 0;
            printf("Fail :: input not match %s %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL ,10, tmp1));
        }
    }

    //return
    if(stat) {
        mpz_set(rop, pf);
    } else {
        mpz_set_ui(rop, 0);
    }
    mpz_clears(inv, tmp1, pf, NULL);
    return stat;
}

//Input: C1x(t) C1y(t) C2x(t) C2y(t); r, z
//Params: r, z
//Output: (C1x+C2x)(t) (C1y+C2y)(t); z
int sum_check_cipher_add(mpz_t rop, mpz_t *V, const mpz_t ri, const mpz_t r, const mpz_t z)
{
    //init
    int stat = 1;
    mpz_t tmp1, poly[2];
    mpz_inits(tmp1, poly[0], poly[1], NULL);

    //Prover
    //update V :: from lsb
    mod_1neg(tmp1, z);
    for(int i = 0; i < 2; i ++){
        mod_mult(V[2 * i], V[2 * i], tmp1);
        mod_mult(V[2 * i + 1], V[2 * i + 1], z);
        mod_add(poly[i], V[2 * i], V[2 * i + 1]);
    }

    //Verifier
    mod_add(tmp1, poly[0], poly[1]);
    if(mpz_cmp(tmp1, ri)) {
        stat = 0;
        printf("Fail :: 1st sum check %s %s\n", mpz_get_str(NULL, 10, tmp1), mpz_get_str(NULL, 10, ri));
    }

    //return
    if(stat) {
        mod_1neg(tmp1, r);
        mod_mult(poly[0], poly[0], tmp1);
        mod_mult(poly[1], poly[1], r);
        mod_add(rop, poly[0], poly[1]);
    } else { 
        mpz_set_ui(rop, 0);
    }
    
    //free
    mpz_clears(tmp1, poly[0], poly[1], NULL);

    return stat;
}

//Input: Cx(t) Cy(t); z
//Params: logN, l1, l2, r[1 + logN + logQ], z[1 + logN + logQ]
//Output: [Cx](t) [Cy](t); z
int sum_check_cipher_rescale(mpz_t rop, mpz_t *C, const mpz_t ri, const int logN, const int l1, const int l2, const mpz_t *rn, const mpz_t *rd, const mpz_t *rb, const mpz_t val, const mpz_t *zn, const mpz_t *zdb)
{
    //ln = 1
    //ld = logN
    //lb = logQ
    return sum_check_poly_rounding(rop, C, ri, 1, logN, l1, l2, rn, rd, rb, val, zn, zdb);
}


