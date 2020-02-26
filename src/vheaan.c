#include "vheaan.h"

int sum_check_verification(mpz_t rop, mpz_t **poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);

    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check\n%s\n%s\n%s\n", str, mpz_get_str(NULL, 10, tmp), mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, PRIME));
        stat = 0;
    }
    
    if(stat) {
        polynomial_extrapolate_N(extrap_val, r[num_rounds - 1], poly[0], degree);
        for(int round = 1; round < num_rounds; round ++) {
            mod_add(tmp, poly[round][0], poly[round][1]);
            if(mpz_cmp(tmp, extrap_val)) {
                printf("Fail sumcheck %s :: %dth round %s %s\n", str, round, mpz_get_str(NULL, 10, tmp), mpz_get_str(NULL, 10, extrap_val));
                stat = 0;
                break;
            } else {
                polynomial_extrapolate_N(extrap_val, r[num_rounds - 1 - round], poly[round], degree);
            }
        }
    }

    if(stat)
        mpz_set(rop, extrap_val);
    else
        mpz_set_ui(rop, ERROR);

    mpz_clears(tmp, extrap_val, NULL);

    return stat;
}

int gkr_cipher_mult(
    const mpz_t *V_0r, const mpz_t *V_0c, const mpz_t *V_0d, 
    const mpz_t *V_1l, const mpz_t *V_1r, const mpz_t *V_1c, const mpz_t * V_1d,
    const mpz_t *V_2, 
    const mpz_t *r_0r, const mpz_t *r_0c, const mpz_t *r_0d, const mpz_t *r_0b,
    const mpz_t *r_1l, const mpz_t *r_1r, const mpz_t *r_1c, const mpz_t *r_1d, const mpz_t *r_1b, 
    const mpz_t *r_2, 
    const mpz_t P, const mpz_t EVK1, const mpz_t EVK2, 
    const int log_num, const int bits,
    const mpz_t val
)
{
    return 0;
}


int gkr_cipher_rescale
(
    const mpz_t *V_r_in, const mpz_t *V_d_in,
    const mpz_t *V_c_in,
    const mpz_t *r_c, const mpz_t *r_b,
    const int log_num, const int bef_bits, const int aft_bits,
    const mpz_t val
)
{
    digit_rep = 16;
    int log_bits = 0;
    while ((bef_bits - 1) >> ++ log_bits);
    int dif_bits = bef_bits - aft_bits, ubits = 1 << log_bits, num = 1 << log_num, stat = 1;
    uint64 num_terms;

    mpz_t tmp1, tmp2, scale, cr, a1r, a2r, bhr, bbr, tr, vdr, vrr;
    mpz_inits(tmp1, tmp2, scale, cr, a1r, a2r, bhr, bbr, tr, vdr, vrr, 0);

    mpz_t bh[3], bb[4], c[4], a1[3], a2[3], t[3];
    for (int i = 0; i < 4; i ++)
        mpz_inits(bb[i], c[i], 0);
    for (int i = 0; i < 3; i ++)
        mpz_inits(bh[i], a1[i], a2[i], t[i], 0);
    
    mpz_t *beta_h = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *beta_b = (mpz_t *) malloc(sizeof(mpz_t) * N * num * ubits * 2),
          *alpha1 = (mpz_t *) malloc(sizeof(mpz_t) * ubits),
          *alpha2 = (mpz_t *) malloc(sizeof(mpz_t) * ubits),
          *tau = (mpz_t *) malloc(sizeof(mpz_t) * N), 
          *V_c = (mpz_t *) malloc(sizeof(mpz_t) * N * num * ubits * 2),
          **poly_b = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
          **poly_d = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
          **poly_r = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
          *tmps1 = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num + logN + log_bits)),
          *tmps2 = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num + logN + log_bits)),
          *r_d = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num));

    for (int i = 0; i < (int) num * 2; i ++)
        mpz_inits(beta_h[i], 0);
    for (int i = 0; i < (int) N * num * ubits * 2; i ++)
        mpz_inits(beta_b[i], V_c[i], 0);
    for (int i = 0; i < ubits; i ++)
        mpz_inits(alpha1[i], alpha2[i], 0);
    for(int i = 0; i < (int) N; i ++)
        mpz_init(tau[i]);
    for (int i = 0; i < (int) (1 + log_num + logN + log_bits); i ++) {
        poly_b[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
        poly_r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
        for (int j = 0; j < 3; j ++) {
            mpz_init_set_ui(poly_d[i][j], 0);
            mpz_init_set_ui(poly_r[i][j], 0);
        }
        for (int j = 0; j < 4; j ++)
            mpz_init_set_ui(poly_b[i][j], 0);
    }
    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++)
        mpz_inits(tmps1[i], tmps2[i], 0);
    for (int i = 0; i < 1 + log_num; i ++) {
        mpz_set(tmps2[i], r_b[i + logN + log_bits]);
    }
    printf("GKR Rescale: Memory allocated\n");


    //init mle
    initialize_beta(beta_b, log_num + 1 + logN + log_bits, r_b);
    initialize_beta(beta_h, log_num + 1, tmps2);
    initialize_alpha(alpha1, log_bits, bef_bits);
    initialize_alpha(alpha2, log_bits, dif_bits);
    initialize_tau(tau, logN, val);
    for (int i = 0; i < (int) N * num * 2 * (1 << log_bits); i ++)
        mpz_init_set(V_c[i], V_c_in[i]);
    mod_inv_pow2(scale, dif_bits);
    printf("GKR Rescale: MLE initialized\n");
 
 
    //eval mle
    evaluate_V(vrr, V_r_in, 1 + log_num, tmps2);
    evaluate_V(vdr, V_d_in, 1 + log_num, tmps2);
    evaluate_V(cr, V_c, 1 + log_num + logN + log_bits, r_c);
    evaluate_beta(bbr, r_b, r_c, log_num + 1 + logN + log_bits);

    for (int i = 0; i < 1 + log_num; i ++)
        mpz_set(tmps1[i], r_c[i + log_bits + logN]);
    evaluate_beta(bhr, tmps2, tmps1, log_num + 1);
    
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(tmps1[i], r_c[i + log_bits]);
    evaluate_tau(tr, val, tmps1, logN);
    
    for (int i = 0; i < log_bits; i ++)
        mpz_set(tmps1[i], r_c[i]);
    evaluate_alpha(a1r, bef_bits, tmps1, log_bits);
    evaluate_alpha(a2r, dif_bits, tmps1, log_bits);
    printf("GKR Rescale: MLE Precomputed\n");

   
    //Prove
    num_terms = 1ULL << (1 + log_num + logN + log_bits);
    for (int round = 0; round < (int) (1 + log_num + logN + log_bits); round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++) {

            //linear computation
            mlmap_evaluation_N(bb, 4, num_terms, j, beta_b);
            mlmap_evaluation_N(c, 4, num_terms, j, V_c);

            mlmap_evaluation_N(bh, 3, num_terms >> (logN + log_bits), j >> (logN + log_bits), beta_h);
            mlmap_evaluation_N(t, 3, (num_terms >> log_bits) & (2 * N - 1), (j >> log_bits) & (N - 1), tau);
            mlmap_evaluation_N(a1, 3, num_terms & (2 * ubits - 1), j & (ubits - 1), alpha1);
            mlmap_evaluation_N(a2, 3, num_terms & (2 * ubits - 1), j & (ubits - 1), alpha2);

            //update poly_b
            for (int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, bb[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }

            //update poly_d
            for (int k = 0; k < 3; k ++) {
                mod_mult(tmp1, bh[k], t[k]);
                mod_mult(tmp1, tmp1, a1[k]);
                mod_mult(tmp1, tmp1, c[k]);
                mod_add(poly_d[round][k], poly_d[round][k], tmp1);
            }

            //update poly_r
            for (int k = 0; k < 3; k ++) {
                mod_sub(tmp1, a1[k], a2[k]);
                mod_mult(tmp1, tmp1, bh[k]);
                mod_mult(tmp1, tmp1, t[k]);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, scale);
                mod_add(poly_r[round][k], poly_r[round][k], tmp1);
            }
        }

        //update mle
        num_terms >>= 1;
        
        update_V(V_c, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
        update_V(beta_b, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
        if(round < 1 + log_num) {
            update_V(beta_h, num_terms >> (logN + log_bits), r_c[1 + log_num + logN + log_bits - round - 1]);
        } else if (round < (int) (1 + log_num + logN)) {
            update_V(tau, num_terms >> log_bits, r_c[1 + log_num + logN + log_bits - round - 1]);
        } else {
            update_V(alpha1, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
            update_V(alpha2, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
        }
    }
    printf("GKR Prove: Proof generated\n");

    //Verify
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, 1 + log_num + logN + log_bits, r_c, "Vb") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, bbr);
    if(stat && mpz_cmp(tmp1, tmp2)) { 
        stat = 0;
        printf("Fail :: Vb does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
    }
    if(stat)
        printf("Vb Verified\n");

    stat = stat ? sum_check_verification(tmp2, poly_r, vrr, 3, 1 + log_num + logN + log_bits, r_c, "Vr") : 0;
    mod_sub(tmp1, a1r, a2r);
    mod_mult(tmp1, tmp1, bhr);
    mod_mult(tmp1, tmp1, tr);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, scale);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Vr does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
    }
    if(stat)
        printf("Vr Verified\n");

    stat = stat ? sum_check_verification(tmp2, poly_d, vdr, 3, 1 + log_num + logN + log_bits, r_c, "Vd") : 0;
    mod_mult(tmp1, tr, bhr);
    mod_mult(tmp1, tmp1, a1r);
    mod_mult(tmp1, tmp1, cr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: Vd does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
    }
    if(stat)
        printf("Vd Verified\n");
    printf("GKR Verify: Verification complete (%d)\n", stat);

    return stat;
}


//int gkr_cipher_mult
//(
//    const mpz_t ri_d, const mpz_t ri_r,
//    const mpz_t *V_c_in, 
//    const mpz_t *r_c, const mpz_t *r_b,
//    const int log_num, const int log_bits, const int bef_bits, const int aft_bits,
//    const mpz_t val
//)
//{
//    digit_rep = 16;
//    int dif_bits = bef_bits - aft_bits, ubits = 1 << log_bits, num = 1 << log_num, stat = 1;
//    uint64 num_terms;
//
//    mpz_t tmp1, tmp2, scale, cr, a1r, a2r, bdr, brr, bbr, tr;
//    mpz_inits(tmp1, tmp2, scale, cr, a1r, a2r, bdr, brr, bbr, tr, 0);
//
//    mpz_t bd[3], br[3], bb[4], c[4], a1[3], a2[3], t[3];
//    for (int i = 0; i < 4; i ++)
//        mpz_inits(bb[i], c[i], 0);
//    for (int i = 0; i < 3; i ++)
//        mpz_inits(bd[i], br[i], a1[i], a2[i], t[i], 0);
//    
//    mpz_t **beta = (mpz_t *) malloc(sizeof(mpz_t) * N * num * ubits * 2),
//          *alpha1 = (mpz_t *) malloc(sizeof(mpz_t) * ubits),
//          *alpha2 = (mpz_t *) malloc(sizeof(mpz_t) * ubits),
//          *tau = (mpz_t *) malloc(sizeof(mpz_t) * N), 
//          *V_c = (mpz_t *) malloc(sizeof(mpz_t) * N * num * ubits * 2),
//          **poly_b = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
//          **poly_d = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
//          **poly_r = (mpz_t **) malloc(sizeof(mpz_t *) * (1 + log_num + logN + log_bits)),
//          *tmps = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num + logN + log_bits));
//
//    for (int i = 0; i < (int) N * num * ubits * 2; i ++)
//        mpz_inits(beta[i], V_c[i], 0);
//    for (int i = 0; i < ubits; i ++)
//        mpz_inits(alpha1[i], alpha2[i], 0);
//    for(int i = 0; i < (int) N; i ++)
//        mpz_init(tau[i]);
//    for (int i = 0; i < (int) (1 + log_num + logN + log_bits); i ++) {
//        poly_b[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
//        poly_d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        for (int j = 0; j < 3; j ++) {
//            mpz_init_set_ui(poly_d[i][j], 0);
//            mpz_init_set_ui(poly_r[i][j], 0);
//        }
//        for (int j = 0; j < 4; j ++)
//            mpz_init_set_ui(poly_b[i][j], 0);
//    }
//    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++)
//        mpz_init(tmps[i]);
//    printf("GKR Rescale: Memory allocated\n");
//
//
//    //init mle
//    initialize_beta(beta, log_num + 1 + logN + log_bits, r_b);
//    initialize_alpha(alpha1, log_bits, bef_bits);
//    initialize_alpha(alpha2, log_bits, dif_bits);
//    initialize_tau(tau, logN, val);
//    for (int i = 0; i < (int) N * num * 2 * (1 << log_bits); i ++)
//        mpz_init_set(V_c[i], V_c_in[i]);
//    mod_inv_pow2(scale, dif_bits);
//    printf("GKR Rescale: MLE initialized\n");
// 
// 
//    //eval mle
//    evaluate_V(cr, V_c, 1 + log_num + logN + log_bits, r_c);
//    evaluate_beta(bbr, r_b, r_c, log_num + 1 + logN + log_bits);
//
//    for (int i = 0; i < 1 + log_num; i ++)
//        mpz_set(tmps[i], r_c[i + log_bits + logN]);
//    for (int i = 0; i < 
//    evaluate_beta(bhr, r_d, tmps, log_num + 1);
//    
//    for (int i = 0; i < (int) logN; i ++)
//        mpz_set(tmps[i], r_c[i + log_bits]);
//    evaluate_tau(tr, val, tmps, logN);
//    
//    for (int i = 0; i < log_bits; i ++)
//        mpz_set(tmps[i], r_c[i]);
//    evaluate_alpha(a1r, bef_bits, tmps, log_bits);
//    evaluate_alpha(a2r, dif_bits, tmps, log_bits);
//    printf("GKR Rescale: MLE Precomputed\n");
//
//   
//    //Prove
//    num_terms = 1ULL << (1 + log_num + logN + log_bits);
//    for (int round = 0; round < (int) (1 + log_num + logN + log_bits); round ++) {
//        for(uint64 j = 0; j < num_terms / 2; j ++) {
//
//            //linear computation
//            mlmap_evaluation_N(bb, 4, num_terms, j, beta_b);
//            mlmap_evaluation_N(c, 4, num_terms, j, V_c);
//
//            mlmap_evaluation_N(bd, 3, num_terms >> (logN + log_bits), j >> (logN + log_bits), beta_d);
//            mlmap_evaluation_N(br, 3, num_terms >> (logN + log_bits), j >> (logN + log_bits), beta_r);
//            mlmap_evaluation_N(t, 3, (num_terms >> log_bits) & (2 * N - 1), (j >> log_bits) & (N - 1), tau);
//            mlmap_evaluation_N(a1, 3, num_terms & (2 * ubits - 1), j & (ubits - 1), alpha1);
//            mlmap_evaluation_N(a2, 3, num_terms & (2 * ubits - 1), j & (ubits - 1), alpha2);
//
//            //update poly_b
//            for (int k = 0; k < 4; k ++) {
//                mpz_sub_ui(tmp1, c[k], 1);
//                mod_mult(tmp1, tmp1, c[k]);
//                mod_mult(tmp1, tmp1, bb[k]);
//                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
//            }
//
//            //update poly_d
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bd[k], t[k]);
//                mod_mult(tmp1, tmp1, a1[k]);
//                mod_mult(tmp1, tmp1, c[k]);
//                mod_add(poly_d[round][k], poly_d[round][k], tmp1);
//            }
//
//            //update poly_r
//            for (int k = 0; k < 3; k ++) {
//                mod_sub(tmp1, a1[k], a2[k]);
//                mod_mult(tmp1, tmp1, br[k]);
//                mod_mult(tmp1, tmp1, t[k]);
//                mod_mult(tmp1, tmp1, c[k]);
//                mod_mult(tmp1, tmp1, scale);
//                mod_add(poly_r[round][k], poly_r[round][k], tmp1);
//            }
//        }
//
//        //update mle
//        num_terms >>= 1;
//        
//        update_V(V_c, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
//        update_V(beta_b, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
//
//        if(round < 1 + log_num) {
//            update_V(beta_d, num_terms >> (logN + log_bits), r_c[1 + log_num + logN + log_bits - round - 1]);
//            update_V(beta_r, num_terms >> (logN + log_bits), r_c[1 + log_num + logN + log_bits - round - 1]);
//        } else if (round < (int) (1 + log_num + logN)) {
//            update_V(tau, num_terms >> log_bits, r_c[1 + log_num + logN + log_bits - round - 1]);
//        } else {
//            update_V(alpha1, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
//            update_V(alpha2, num_terms, r_c[1 + log_num + logN + log_bits - round - 1]);
//        }
//    }
//    printf("GKR Prove: Proof generated\n");
//
//    //Verify
//    mpz_set_ui(tmp1, 0);
//    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, 1 + log_num + logN + log_bits, r_c, "Vb") : 0;
//    mpz_sub_ui(tmp1, cr, 1);
//    mod_mult(tmp1, tmp1, cr);
//    mod_mult(tmp1, tmp1, bbr);
//    if(stat && mpz_cmp(tmp1, tmp2)) { 
//        stat = 0;
//        printf("Fail :: Vb does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
//    }
//    if(stat)
//        printf("Vb Verified\n");
//
//    stat = stat ? sum_check_verification(tmp2, poly_r, ri_r, 3, 1 + log_num + logN + log_bits, r_c, "Vr") : 0;
//    mod_sub(tmp1, a1r, a2r);
//    mod_mult(tmp1, tmp1, brr);
//    mod_mult(tmp1, tmp1, tr);
//    mod_mult(tmp1, tmp1, cr);
//    mod_mult(tmp1, tmp1, scale);
//    if(stat && mpz_cmp(tmp1, tmp2)) {
//        stat = 0;
//        printf("Fail :: Vr does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
//    }
//    if(stat)
//        printf("Vr Verified\n");
//
//    stat = stat ? sum_check_verification(tmp2, poly_d, ri_d, 3, 1 + log_num + logN + log_bits, r_c, "Vd") : 0;
//    mod_mult(tmp1, tr, bdr);
//    mod_mult(tmp1, tmp1, a1r);
//    mod_mult(tmp1, tmp1, cr);
//    if(stat && mpz_cmp(tmp1, tmp2)) {
//        stat = 0;
//        printf("Fail :: Vd does not match %s %s\n", mpz_get_str(NULL, digit_rep, tmp1), mpz_get_str(NULL, digit_rep, tmp2));
//    }
//    if(stat)
//        printf("Vd Verified\n");
//    printf("GKR Verify: Verification complete (%d)\n", stat);
//
//
//    //free
//    mpz_clears(tmp1, tmp2, scale, cr, a1r, a2r, bdr, brr, bbr, tr, 0);
//    for (int i = 0; i < 4; i ++)
//        mpz_clears(bb[i], c[i], 0);
//    for (int i = 0; i < 3; i ++)
//        mpz_clears(bd[i], br[i], a1[i], a2[i], t[i], 0);
//    for (int i = 0; i < (int) num * 2; i ++)
//        mpz_clears(beta_d[i], beta_r[i], 0);
//    for (int i = 0; i < (int) N * num * ubits * 2; i ++)
//        mpz_clears(beta_b[i], V_c[i], 0);
//    for (int i = 0; i < ubits; i ++)
//        mpz_clears(alpha1[i], alpha2[i], 0);
//    for(int i = 0; i < (int) N; i ++)
//        mpz_clear(tau[i]);
//    for (int i = 0; i < (int) (1 + log_num + logN + log_bits); i ++) {
//        for (int j = 0; j < 3; j ++)
//            mpz_clears(poly_d[i][j], poly_r[i][j], 0);
//        for (int j = 0; j < 4; j ++)
//            mpz_clear(poly_b[i][j]);
//        free(poly_b[i]);
//        free(poly_d[i]);
//        free(poly_r[i]);
//    }
//    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits); i ++)
//        mpz_clear(tmps[i]);
//    free(tmps);
//    free(poly_b);
//    free(poly_d);
//    free(poly_r);
//    free(V_c);
//    free(tau);
//    free(alpha1);
//    free(alpha2);
//    free(beta_d);
//    free(beta_r);
//    free(beta_b);
//    printf("GKR free: Memory free\n");
//    return stat;
//}


//int gkr_cipher_mult()
//{};

//int gkr_cipher_add()
//{};

//int gkr_cipher_cmult()
//{};

/*
int gkr_cipher_rescale
(
    mpz_t rop,
    mpz_t *C, mpz_t *Vd, mpz_t *Vr,
    const mpz_t *r,
    const mpz_t *z, const mpz_t *z_bits,
    const int log_num, const int log_bits, const int bef_bits, const int aft_bits,
    const mpz_t *commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t *sks
)
{
    clock_t t clock();
    int stat = 1;

    //ephemeral
    mpz_t tmp1, tmp2, cr, a1r, a2r, br, b[4], c[4], a1[3], a2[3];
    mpz_t *z_all;

    //commit
    evaluate_V(cr, c, log_num + log_bef_bits, r);

    //alpha
    mpz_t *r_bits = (mpz_t *) malloc(log_bits * sizeof(mpz_t)),
          *r_num = (mpz_t *) malloc((log_num + 1) * sizeof(mpz_t));
    evaluate_alpha(a1r, bef_bits, r_bits, log_bits);
    evaluate_alpha(a2r, bef_bits - aft_bits, r_bits, log_bits);
    mod_inv_pow2(tmp1, bef-bits - aft_bits);
    mod_sub(a2r, a1r, a2r);
    mod_mult(a2r, a2r, tmp1);

    //beta
    evaluate_beta(br, z, r, log_bits + log_num + 1);

    //mlmap
    mpz_t *alpha1 = (mpz_t *) malloc(sizeof(mpz_t) << log_bits),
          *alpha2 = (mpz_t *) malloc(sizeof(mpz_t) << log_bits),
          *beta = (mpz_t *) malloc(sizeof(mpz_t) << (log_bits + log_num + 1));
    initialize_alpha(alpha1, log_bits, bef_bits);
    initialize_alpha(alpha2, log_bits, aft_bits);
    initialize_beta(beta, z_all, r, log_num + 1);

    //polynomial
    mpz_t **poly_b = (mpz_t **) malloc(sizeof(mpz_t *) * (log_bits + log_num + 1)),
          **poly_d = (mpz_t **) malloc(sizeof(mpz_t *) * log_bits);
          **poly_r = (mpz_t **) malloc(sizeof(mpz_t *) * log_bits);

    //prove
    uint64 num_terms = 1ULL << (log_bits + log_num + 1);
    for (int round = 0; round < log_bits + log_num + 1; round ++) {
        for (uint64 j = 0; j < num_terms / 2; j ++) {
            
            //lin eval
            if (round < log_bits) {
            } else {

            }

            mlmap_evaluation_N(b, 4, num_terms, j, beta);
            mlmap_evaluation_N(c, 4, num_terms, j, C);

            //upd poly
            for (int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }
            for (int k = 0; k < )
                poly_d
            for (int k)
                poly_r
        }

        num_terms >>= 1;

        //update V
        if (round < log_bits) {
            update_V (beta1)
        } else {

        }

    }

}

int sum_check_cipher_mult
(
    mpz_t rop, mpz_t *CR1_in, mpz_t *CR2_in, mpz_t *V3, mpz_t *V2, mpz_t *V1l, mpz_t *V1rr, mpz_t *VR1x, mpz_t *VR1y, mpz_t *VR2x, mpz_t *VR2y,
    const mpz_t ri_in,
    const mpz_t P, const mpz_t *evk, const mpz_t val,
    const int logN, const int bit,
    const mpz_t r, const mpz_t *z,
    const mpz_t *zR1d, const mpz_t *zR1b, const mpz_t *rR1d, const mpz_t *rR1b,
    const mpz_t *zR2d, const mpz_t *zR2b, const mpz_t *rR2d, const mpz_t *rR2b,
	const mpz_t z1,
	const mpz_t* commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t* sks,
	const mpz_t* commits2, const uint64 lck2, const uint64 key_num2, const uint64 test_comm_num2, const mpz_t* sks2, const mpz_t gen
)
{
    //Init
    double t_V = 0, t_P = 0, t_dummy = 0, time_diff;
    int lb = 0;
    while((bit - 1) >> ++ lb);
    int stat = 1;
	clock_t timer = clock();
    uint64 num_terms, N = 1ULL << logN, B = 1ULL << lb;
    mpz_t tmp1, tmp2, tN, ri, cr, br, a1r, a2r, a3r, tr, g[4], a1[4], a2[4], a3[4], b[4], t[4], c[4], v[4], v1[4], v2[4], ri0, ri1, rir, ril, u[2];
    mpz_t **poly_b, **poly_d1, **poly_d2, **poly_r, *R, *Z, *RR;
    mpz_t *alpha1, *alpha2, *alpha3, *tau, *beta, *gamma, *CR1, *CR2;
    
    mpz_inits(u[0], u[1], NULL);
    mpz_inits(tmp1, tmp2, tN, ri, cr, br, a1r, a2r, a3r, tr, ri0, ri1, rir, ril, NULL);
    mpz_set(ri, ri_in);

    poly_b = (mpz_t**) malloc(sizeof(mpz_t*) * (1 + logN + lb + 2));
    poly_d1 = (mpz_t**) malloc(sizeof(mpz_t*) * logN);
    poly_d2 = (mpz_t**) malloc(sizeof(mpz_t*) * (lb + 2));
    poly_r = (mpz_t**) malloc(sizeof(mpz_t*) * (logN + lb + 2));

    CR1 = (mpz_t*) malloc(sizeof(mpz_t) << (1 + logN + lb + 2));
    CR2 = (mpz_t*) malloc(sizeof(mpz_t) << (logN + lb + 2));
    
    R = (mpz_t*) malloc(sizeof(mpz_t) * (1 + logN + lb + 2));
	RR = (mpz_t*)malloc(sizeof(mpz_t) * (1 + logN + lb + 2));
    Z = (mpz_t*) malloc(sizeof(mpz_t) * (1 + logN + lb + 2));
    alpha1 = (mpz_t*) malloc(sizeof(mpz_t) << (lb + 2));
    alpha2 = (mpz_t*) malloc(sizeof(mpz_t) << (lb + 2));
    alpha3 = (mpz_t*) malloc(sizeof(mpz_t) << (lb + 2));
    tau = (mpz_t*) malloc(sizeof(mpz_t) << logN);
    beta = (mpz_t*) malloc(sizeof(mpz_t) << (1 + logN + lb + 2));
    gamma = (mpz_t*) malloc(sizeof(mpz_t) * 4);

	for (int i = 0; i < 1 + logN + lb + 2; i++)
		mpz_init(RR[i]);

    for(int i = 0; i < 4; i ++)
        mpz_inits(g[i], a1[i], a2[i], a3[i], b[i], t[i], c[i], v[i], v1[i], v2[i], NULL);
    
    for(uint64 i = 0; i < (1ULL << (1 + logN + lb + 2)); i ++){
        mpz_init_set(CR1[i], CR1_in[i]);
    }

    for(uint64 i = 0; i < (1ULL << (logN + lb + 2)); i ++){
        mpz_init_set(CR2[i], CR2_in[i]);
    }

    for(int i = 0; i < logN; i ++) {
        poly_d1[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        for(int j = 0; j < 4; j ++)
            mpz_init(poly_d1[i][j]);
    }
    for(int i = 0; i < lb + 2; i ++) {
        poly_d2[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        for(int j = 0; j < 4; j ++)
            mpz_init(poly_d2[i][j]);
    }
    for(int i = 0; i < 1 + logN + lb + 2; i ++) {
        poly_b[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        for(int j = 0; j < 4; j ++)
            mpz_init(poly_b[i][j]);
        mpz_inits(R[i], Z[i], NULL);
    }
    for(int i = 0; i < logN + lb + 2; i ++) {
        poly_r[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        for(int j = 0; j < 4; j ++)
            mpz_init(poly_r[i][j]);
    }

    for(uint64 i = 0; i < (1ULL << (lb + 2)); i ++)
        mpz_inits(alpha1[i], alpha2[i], alpha3[i], NULL);
    for(uint64 i = 0; i < (1ULL << logN); i ++)
        mpz_init(tau[i]);
    for(uint64 i = 0; i < (1ULL << (1 + logN + lb + 2)); i ++)
        mpz_init(beta[i]);
    for(uint64 i = 0; i < 4; i ++)
        mpz_init(gamma[i]);
    mod_pow(tN, val, N);
    printf("Init done\n");

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_dummy += time_diff;
	printf("%f ms dummy for Init \n", time_diff);
	timer = clock();


    //LayerR1

    ////LayerR1 init
    for(int i = 0; i < lb + 2; i ++) {
        mpz_set(R[i], rR1b[i]);
        mpz_set(Z[i], zR1b[i]);
    }
    for(int i = 0; i < logN; i ++) {
        mpz_set(R[i + lb + 2], rR1d[i]);
        mpz_set(Z[i + lb + 2], zR1d[i]);
    }
    for(int i = 0; i < 1; i ++) {
        mpz_set(R[i + logN + lb + 2], z[i]);
        mpz_set(Z[i + logN + lb + 2], z[i]);
    }
    for(int i = 0; i < 1 + logN + lb + 2; i ++)
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_b[i][k], 0);
    for(int i = 0; i < logN + lb + 2; i ++)
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_r[i][k], 0);
    for(int i = 0; i < logN; i ++) 
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_d1[i][k], 0);
    for(int i = 0; i < lb + 2; i ++) 
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_d2[i][k], 0);

    printf("LayerR1 init done\n");
    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_dummy += time_diff;
	printf("%f ms dummy for init \n", time_diff);
	timer = clock();

   
    ////LayerR1 Prover init
    initialize_beta(beta, 1 + logN + (lb + 2), Z);
    initialize_alpha(alpha1, lb + 2, bit * 4);
    initialize_alpha(alpha2, lb + 2, bit * 2);
    initialize_alpha(alpha3, lb + 2, bit * 1);
    initialize_tau(tau, logN, val);
    printf("LayerR1 prover init done\n");

    ////LayerR1 Proof
    num_terms = 1ULL << (1 + logN + (lb + 2));
    for(int round = 0; round < 1 + logN + (lb + 2); round ++) {

        for(uint64 j = 0; j < num_terms / 2; j ++) {

            //linear computation
            mlmap_evaluation_N(c, 4, num_terms, j, CR1);
            mlmap_evaluation_N(b, 4, num_terms, j, beta);

            if(round < 1) {

            } else if (round < 1 + logN) {
                
                mlmap_evaluation_N(t, 3, num_terms >> (lb + 2), j >> (lb + 2), tau);
                if((j & ((1 << (lb +2)) - 1)) == 0){
                    mlmap_evaluation_N(v1, 3, num_terms >> (lb + 2), j >> (lb + 2), VR1x);
                    mlmap_evaluation_N(v2, 3, num_terms >> (lb + 2), j >> (lb + 2), VR1y);
                }
                mlmap_evaluation_N(a2, 3, 0, j & ((1 << (lb + 2)) - 1), alpha2);
                mlmap_evaluation_N(a3, 3, 0, j & ((1 << (lb + 2)) - 1), alpha3);

            } else {

                mlmap_evaluation_N(t, 3, 0, 0, tau);
                mlmap_evaluation_N(a1, 3, num_terms, j, alpha1);
                mlmap_evaluation_N(a2, 3, num_terms, j, alpha2);
                mlmap_evaluation_N(a3, 3, num_terms, j, alpha3);
            }

            //update poly_b
            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }

            //update poly_d
            if(round < 1){
            }else if(round < 1 + logN) {
                if((j & ((1 << (lb +2)) - 1)) == 0){
                    for(int k = 0; k < 3; k ++) {
                        mod_mult(tmp1, tN, v2[k]);
                        mod_add(tmp1, tmp1, v1[k]);
                        mod_mult(tmp1, tmp1, t[k]);
                        mod_add(poly_d1[round - 1][k], poly_d1[round - 1][k], tmp1);
                    }
                } else {
                }
            } else {
                for(int k = 0; k < 3; k ++) {
                    mod_mult(tmp1, c[k], a1[k]);
                    mod_add(poly_d2[round - 1 - logN][k], poly_d2[round - 1 - logN][k], tmp1);
                }
            }

            //update poly_r
            if(round < 1 ) {
            } else {
                for(int k = 0; k < 3; k ++) {
                    mod_sub(tmp1, a2[k], a3[k]);
                    mod_mult(tmp1, tmp1, t[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_r[round - 1][k], poly_r[round - 1][k], tmp1);
                }
            }
        }

        num_terms >>= 1;
        update_V(CR1, num_terms, R[1 + logN + lb + 2 - round - 1]);
        update_V(beta, num_terms, R[1 + logN + lb + 2 - round - 1]);
        if(round < 1) {
            update_V(VR1x, num_terms >> (lb + 2), R[1 + logN + lb + 2 - 1 - round]);
            update_V(VR1y, num_terms >> (lb + 2), R[1 + logN + lb + 2 - 1 - round]);
        } else if (round < 1 + logN) {
            update_V(tau, num_terms >> (lb + 2), R[1 + logN + lb + 2 - 1 - round]);
            update_V(VR1x, num_terms >> (lb + 2), R[1 + logN + lb + 2 - 1 - round]);
            update_V(VR1y, num_terms >> (lb + 2), R[1 + logN + lb + 2 - 1 - round]);
        } else {
            update_V(alpha1, num_terms, R[1 + logN + lb + 2 - 1 - round]);
            update_V(alpha2, num_terms, R[1 + logN + lb + 2 - 1 - round]);
            update_V(alpha3, num_terms, R[1 + logN + lb + 2 - 1 - round]);
        }
    }
    mpz_set(ri0, VR1x[0]);
    mpz_set(ri1, VR1y[0]);

    printf("LayerR1 proof done\n");
    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover prove R1\n", time_diff);
	timer = clock();



    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_dummy += time_diff;
	printf("%f ms dummy for reversing \n", time_diff);
	timer = clock();


    ////LayerR1 commit
    evaluate_beta(br, Z, R, 1 + logN + lb + 2);
    evaluate_alpha(a1r, 4 * bit, rR1b, lb + 2);
    evaluate_alpha(a2r, 2 * bit, rR1b, lb + 2);
    evaluate_alpha(a3r, bit, rR1b, lb + 2);
    evaluate_tau(tr, val, rR1d, logN);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier for eval greeks \n", time_diff);
	timer = clock();



	// COMMIT PART
	for (int i = 0; i < 1 + logN + lb + 2; ++i) {									// rr = reverse order of r // for commit only
		mpz_init(RR[i]);
		mpz_set(RR[i], R[1 + logN + lb + 2  - 1 - i]);
	}

	mpz_t* test_outs = (mpz_t*)malloc((key_num) * sizeof(mpz_t));
	mpz_t* test_evalpts = (mpz_t*)malloc((test_comm_num) * sizeof(mpz_t));	// Induced from log (test_comm_num) = l-lck # of end-variables of r

	for (uint64 i = 0; i < test_comm_num; ++i)								// Initialize Evalpts for Test
		mpz_init2(test_evalpts[i], BITS);

	timer = clock();
	kxi_eval(test_evalpts, test_comm_num, 1 + logN + lb + 2 - lck, &(RR[lck]));

	printf("%f ms Verifier to evaluate before commit \n", (double)(clock() - timer) / CLOCKS_PER_SEC * 1000);

	commit_open_binary(test_outs, CR1_in, test_evalpts, commits, key_num, test_comm_num, sks, gen);

	printf("for Open & Check Commit \n");


	for (long i = 0; i < lck; ++i)					// front of rr = same order as r // for commit only
		mpz_set(RR[i], R[1 + logN + lb + 2 - lck + i]);

	timer = clock();
	evaluate_V(cr, test_outs, lck, RR);
	printf("%f ms Verifier to evaluate with commit \n", (double)(clock() - timer) / CLOCKS_PER_SEC * 1000);
	timer = clock();
    printf("LayerR1 commit done\n");


    ////LayerR1 Verifier
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, 1 + logN + lb + 2, R, "LayerR1 Vb") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, br);
    if(stat && mpz_cmp(tmp1, tmp2)) { 
        stat = 0;
        printf("Fail :: LayerR1 Vb does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR1 Vb Verified\n");

    mpz_set_ui(tmp2, 2);
    mod_pow(tmp2, tmp2, bit);
    mod_mult(tmp1, tmp2, ri);
    stat = stat ? sum_check_verification(tmp2, poly_r, tmp1, 3, logN + lb + 2, R, "LayerR1 Vr") : 0;
    mod_sub(tmp1, a2r, a3r);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR1 Vr does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR1 Vr Verified\n");

    mod_add(tmp1, poly_d1[0][1], poly_d1[0][0]);
    stat = stat ? sum_check_verification(tmp2, poly_d1, tmp1, 3, logN, rR1d, "LayerR1 Vd1") : 0;
    mod_mult(tmp1, ri1, tN);
    mod_add(tmp1, tmp1, ri0);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR1 Vd1 does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR1 Vd1 Verified\n");


    mpz_set_ui(tmp1, 2);
    mod_pow(tmp1, tmp1, 4 * bit - 1);
    mod_sub(tmp1, tmp1, ri1);
    mod_add(tmp1, tmp1, ri0);
    stat = stat ? sum_check_verification(tmp2, poly_d2, tmp1, 3, lb + 2, rR1b, "LayerR1 Vd2") : 0;
    mod_mult(tmp1, a1r, cr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR1 Vd2 does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR1 Vd2 Verified\n");


    ////LayerR1 Return
    if(stat) {
        mod_add(ri, poly_d1[0][1], poly_d1[0][0]);
        printf("LayerR1 Verified\n");
    } else {
        mpz_set_ui(ri, 0);
    }
   
    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier for verification \n", time_diff);
	timer = clock();




    //LayerAdd

    ////LayerAdd Prover
    evaluate_V(ril, V1l, 1, z);
    mpz_set(rir, V1rr[0]);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover for prove add \n", time_diff);
	timer = clock();


    ////LayerAdd Verifier
    evaluate_V(tmp1, evk, 1, z);
    mod_mult(tmp1, tmp1, rir);
    mod_add(tmp1, tmp1, ril);//tmp1=ri0 + ri1*((1-z)*evk1+z*evk2)
    if(stat & mpz_cmp(tmp1, ri)) {
        stat = 0;
        printf("Fail :: LayerAdd doees not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, ri));
    }

    ////LayerAdd Return
    if(stat) {
        printf("LayerAdd Verified\n");
    } else {
        mpz_set_ui(ril, 0);
        mpz_set_ui(rir, 0);
    }

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier for add verification  \n", time_diff);
	timer = clock();



   
    //LayerR2

    ////LayerR2 init
    for(int i = 0; i < lb + 2; i ++) {
        mpz_set(R[i], rR2b[i]);
        mpz_set(Z[i], zR2b[i]);
    }
    for(int i = 0; i < logN; i ++) {
        mpz_set(R[i + lb + 2], rR2d[i]);
        mpz_set(Z[i + lb + 2], zR2d[i]);
    }
    for(int i = 0; i < logN + lb + 2; i ++)
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_b[i][k], 0);
    for(int i = 0; i < logN + lb + 2; i ++)
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_r[i][k], 0);
    for(int i = 0; i < logN; i ++) 
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_d1[i][k], 0);
    for(int i = 0; i < lb + 2; i ++) 
        for(int k = 0; k < 4; k ++)
            mpz_set_ui(poly_d2[i][k], 0);

    printf("LayerR2 init done\n");
    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_dummy += time_diff;
	printf("%f ms dummy init R2 \n", time_diff);
	timer = clock();


    ////LayerR2 Prover init
    initialize_beta(beta, 0 + logN + (lb + 2), Z);
    initialize_alpha(alpha1, lb + 2, bit * 4);
    initialize_alpha(alpha2, lb + 2, bit * 1);
    initialize_tau(tau, logN, val);
    printf("LayerR2i prover init done\n");

    ////LayerR2 Proof
    num_terms = 1ULL << (0 + logN + (lb + 2));
    for(int round = 0; round < 0 + logN + (lb + 2); round ++) {

        for(uint64 j = 0; j < num_terms / 2; j ++) {

            //linear computation
            mlmap_evaluation_N(c, 4, num_terms, j, CR2);
            mlmap_evaluation_N(b, 4, num_terms, j, beta);

            if(round < 0) {

            } else if (round < 0 + logN) {
                
                mlmap_evaluation_N(t, 3, num_terms >> (lb + 2), j >> (lb + 2), tau);
                if((j & ((1 << (lb + 2)) - 1)) == 0){
                    mlmap_evaluation_N(v1, 3, num_terms >> (lb + 2), j >> (lb + 2), VR2x);
                    mlmap_evaluation_N(v2, 3, num_terms >> (lb + 2), j >> (lb + 2), VR2y);
                }
                mlmap_evaluation_N(a2, 3, 0, j & ((1 << (lb + 2)) - 1), alpha2);

            } else {

                mlmap_evaluation_N(t, 3, 0, 0, tau);
                mlmap_evaluation_N(a1, 3, num_terms, j, alpha1);
                mlmap_evaluation_N(a2, 3, num_terms, j, alpha2);
            }

            //update poly_b
            for(int k = 0; k < 4; k ++) {
                mpz_sub_ui(tmp1, c[k], 1);
                mod_mult(tmp1, tmp1, c[k]);
                mod_mult(tmp1, tmp1, b[k]);
                mod_add(poly_b[round][k], poly_b[round][k], tmp1);
            }

            //update poly_d
            if(round < 0){
            }else if(round < 0 + logN) {
                if((j & ((1 << (lb +2)) - 1)) == 0){
                    for(int k = 0; k < 3; k ++) {
                        mod_mult(tmp1, tN, v2[k]);
                        mod_add(tmp1, tmp1, v1[k]);
                        mod_mult(tmp1, tmp1, t[k]);
                        mod_add(poly_d1[round - 0][k], poly_d1[round - 0][k], tmp1);
                    }
                } else {
                }
            } else {
                for(int k = 0; k < 3; k ++) {
                    mod_mult(tmp1, c[k], a1[k]);
                    mod_add(poly_d2[round - 0 - logN][k], poly_d2[round - 0 - logN][k], tmp1);
                }
            }

            //update poly_r
            if(round < 0 ) {
            } else {
                for(int k = 0; k < 3; k ++) {
                    mod_mult(tmp1, a2[k], t[k]);
                    mod_mult(tmp1, tmp1, c[k]);
                    mod_add(poly_r[round - 0][k], poly_r[round - 0][k], tmp1);
                }
            }
        }

        num_terms >>= 1;
        update_V(CR2, num_terms, R[0 + logN + lb + 2 - round - 1]);
        update_V(beta, num_terms, R[0 + logN + lb + 2 - round - 1]);
        if(round < 0) {
            update_V(VR1x, num_terms >> (lb + 2), R[0 + logN + lb + 2 - 1 - round]);
            update_V(VR1y, num_terms >> (lb + 2), R[0 + logN + lb + 2 - 1 - round]);
        } else if (round < 0 + logN) {
            update_V(tau, num_terms >> (lb + 2), R[0 + logN + lb + 2 - 1 - round]);
            update_V(VR2x, num_terms >> (lb + 2), R[0 + logN + lb + 2 - 1 - round]);
            update_V(VR2y, num_terms >> (lb + 2), R[0 + logN + lb + 2 - 1 - round]);
        } else {
            update_V(alpha1, num_terms, R[0 + logN + lb + 2 - 1 - round]);
            update_V(alpha2, num_terms, R[0 + logN + lb + 2 - 1 - round]);
        }
    }
    mpz_set(ri0, VR2x[0]);
    mpz_set(ri1, VR2y[0]);
    printf("LayerR2 proof done\n");

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover proves R2 \n", time_diff);
	timer = clock();


    ////LayerR2 commit
    evaluate_beta(br, Z, R, 0 + logN + lb + 2);
    evaluate_alpha(a1r, 4 * bit, rR2b, lb + 2);
    evaluate_alpha(a2r, 1 * bit, rR2b, lb + 2);
    evaluate_tau(tr, val, rR2d, logN);
    printf("LayerR2 commit done\n");

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier eval Greeks \n", time_diff);
	timer = clock();



	// COMMIT PART

	for (int i = 0; i < logN + lb + 2; ++i) {									// rr = reverse order of r // for commit only
		mpz_set(RR[i], R[logN + lb + 2  - 1 - i]);
	}

	mpz_t* test_outs2 = (mpz_t*)malloc((key_num2) * sizeof(mpz_t));
	mpz_t* test_evalpts2 = (mpz_t*)malloc((test_comm_num2) * sizeof(mpz_t));	// Induced from log (test_comm_num) = l-lck # of end-variables of r

	for (uint64 i = 0; i < test_comm_num2; ++i)								// Initialize Evalpts for Test
		mpz_init2(test_evalpts2[i], BITS);

	timer = clock();

	kxi_eval(test_evalpts2, test_comm_num2, logN + lb + 2 - lck2, &(RR[lck2]));

	printf("%f ms Verifier to evaluate before commit \n", (double)(clock() - timer) / CLOCKS_PER_SEC * 1000);

	commit_open_binary(test_outs2, CR2_in, test_evalpts2, commits2, key_num2, test_comm_num2, sks2, gen);

	printf("for Open & Check Commit \n");


	for (long i = 0; i < lck2; ++i)					// front of rr = same order as r // for commit only
		mpz_set(RR[i], R[logN + lb + 2 - lck2 + i]);

	timer = clock();
	evaluate_V(cr, test_outs2, lck2, RR);
	printf("%f ms Verifier to evaluate with commit \n", (double)(clock() - timer) / CLOCKS_PER_SEC * 1000);


	timer = clock();


	printf("LayerR2 commit done\n");


    ////LayerR2 Verifier
    mpz_set_ui(tmp1, 0);
    stat = stat ? sum_check_verification(tmp2, poly_b, tmp1, 4, 0 + logN + lb + 2, R, "LayerR2 Vb") : 0;
    mpz_sub_ui(tmp1, cr, 1);
    mod_mult(tmp1, tmp1, cr);
    mod_mult(tmp1, tmp1, br);
    if(stat && mpz_cmp(tmp1, tmp2)) { 
        stat = 0;
        printf("Fail :: LayerR2 Vb does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR2 Vb Verified\n");

    mpz_set(tmp1, rir);
    stat = stat ? sum_check_verification(tmp2, poly_r, tmp1, 3, logN + lb + 2, R, "LayerR2 Vr") : 0;
    mod_mult(tmp1, a2r, cr);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR2 Vr does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR2 Vr Verified\n");

    mod_add(tmp1, poly_d1[0][1], poly_d1[0][0]);
    stat = stat ? sum_check_verification(tmp2, poly_d1, tmp1, 3, logN, rR2d, "LayerR2 Vd1") : 0;
    mod_mult(tmp1, ri1, tN);
    mod_add(tmp1, tmp1, ri0);
    mod_mult(tmp1, tmp1, tr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR2 Vd1 does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR2 Vd1 Verified\n");

    mpz_set_ui(tmp1, 2);
    mod_pow(tmp1, tmp1, 4 * bit - 1);
    mod_sub(tmp1, tmp1, ri1);
    mod_add(tmp1, tmp1, ri0);
    stat = stat ? sum_check_verification(tmp2, poly_d2, tmp1, 3, lb + 2, rR2b, "LayerR2 Vd2") : 0;
    mod_mult(tmp1, a1r, cr);
    if(stat && mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail :: LayerR2 Vd2 does not match %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, tmp2));
    }
    if(stat)
        printf("LayerR2 Vd2 Verified\n");


    ////LayerR2 Return
    if(stat) {
        mod_add(rir, poly_d1[0][1], poly_d1[0][0]);
        printf("LayerR2 Verified\n");
    } else {
        mpz_set_ui(rir, 0);
    }

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier verifies R2 \n", time_diff);
	timer = clock();

    

    //LayerSplit
    ////LayerSplit Provoer
    mpz_set_ui(Z[0], 0);
    mpz_set(Z[1], z[0]);
    evaluate_V(poly_r[0][0], V2, 2, Z);

    mpz_set_ui(Z[0], 1);
    mod_1neg(Z[1], z[0]);
    evaluate_V(poly_r[0][1], V2, 2, Z);

    mod_1neg(tmp1, z[0]);
    mpz_mul_ui(tmp1, tmp1, 3);
    mpz_sub_ui(tmp1, tmp1, 1);
    mod(Z[1], tmp1);
    mpz_set_ui(Z[0], 2);
    evaluate_V(poly_r[0][2], V2, 2, Z);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover proves Split \n", time_diff);
	timer = clock();

    ////LayerSplit Verifier
    mod_inv(tmp1, P);
    mod_mult(tmp1, tmp1, ril);
    mod_1neg(tmp2, z[0]);
    mod_mult(tmp2, tmp2, rir);

    mod_add(tmp1, tmp1, tmp2);//tmp1<-V1l(z)/P + (1-z)V1r
    mod_add(tmp2, poly_r[0][0], poly_r[0][1]);//tmp2<-f(0)+f(1)
    if(stat & mpz_cmp(tmp1, tmp2)) {
        stat = 0;
        printf("Fail LayerSplit: %s %s\n", mpz_get_str(0, 16, tmp1), mpz_get_str(0, 16, tmp2));
    }

    if(stat) {
        extrapolate(ri, poly_r[0], 3, z1);
        printf("LayerSplit Verified\n");
    } else {
        mpz_set_ui(ri, 0);
    }

    mpz_set(u[0], z1);
    mod_1neg(tmp1, z1);
    mod_1neg(tmp2, z[0]);
    mod_mult(tmp1, tmp1, z[0]);
    mod_mult(tmp2, tmp2, z1);
    mod_add(u[1], tmp1, tmp2);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier verifies Splits \n", time_diff);
	timer = clock();



    //LayerProd
    ////LayerProd Prover
    mpz_set(Z[0], u[1]);
    mpz_set_ui(Z[1], 0);
    evaluate_V(poly_r[0][0], V3, 2, Z);
    
    mpz_set(Z[0], u[0]);
    mpz_set_ui(Z[1], 1);
    evaluate_V(poly_r[0][1], V3, 2, Z);

    mod_add(tmp1, u[0], u[0]);
    mod_sub(Z[0], tmp1, u[1]);
    mpz_set_ui(Z[1], 2);
    evaluate_V(poly_r[0][2], V3, 2, Z);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Proves Prod \n", time_diff);
	timer = clock();

    ////LayerProd Verifier
    mod_mult(tmp1, poly_r[0][0], poly_r[0][1]);
    if(stat & mpz_cmp(tmp1, ri)) {
        stat = 0;
        printf("Fail :: LayerProd %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, ri));
    }

    ////LayerProd Return
    if(stat) {
        printf("LayerProd Verified\n");
        extrapolate(ri, poly_r[0], 3, r);
    } else {
        mpz_set_ui(ri, 0);
    }

    mpz_set(R[1], r);
    mod_1neg(tmp1, r);
    mod_mult(tmp1, tmp1, u[1]);
    mod_mult(tmp2, r, u[0]);
    mod_add(R[0], tmp1, tmp2);
    evaluate_V(tmp1, V3, 2, R);

    
    //Return
    if(stat) {
        printf("Circuit Verified\n");
        mpz_set(rop, ri);
    } else {
        mpz_set_ui(rop, 0);
    }

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier verifies Prod \n", time_diff);
	timer = clock();

    
    
    //Free
    mpz_clears(tmp1, tmp2, tN, ri, cr, br, a1r, a2r, a3r, tr, ri0, ri1, rir, ril, NULL);
    for(int i = 0; i < logN; i ++) {
        for(int j = 0; j < 4; j ++)
            mpz_clear(poly_d1[i][j]);
        free(poly_d1[i]);
    }
    for(int i = 0; i < lb + 2; i ++) {
        for(int j = 0; j < 4; j ++)
            mpz_clear(poly_d2[i][j]);
        free(poly_d2[i]);
    }
    for(int i = 0; i < 1 + logN + lb + 2; i ++) {
        for(int j = 0; j < 4; j ++) {
            mpz_clears(poly_b[i][j], NULL);
        }
        mpz_inits(R[i], Z[i], NULL);
        free(poly_b[i]);
    }
    for(int i = 0; i < logN + lb + 2; i ++) {
        for(int j = 0; j < 4; j ++) {
            mpz_clears(poly_r[i][j], NULL);
        }
        free(poly_r[i]);
    }
    for(uint64 i = 0; i < (1ULL << (1 + logN + lb + 2)); i ++){
        mpz_clear(CR1[i]);
    }

    for(uint64 i = 0; i < (1ULL << (logN + lb + 2)); i ++){
        mpz_clear(CR2[i]);
    }
    free(CR1);
    free(CR2);
    mpz_clears(u[0], u[1], NULL);

    for(uint64 i = 0; i < (1ULL << (lb + 2)); i ++)
        mpz_clears(alpha1[i], alpha2[i], alpha3[i], NULL);
    for(uint64 i = 0; i < (1ULL << logN); i ++)
        mpz_clear(tau[i]);
    for(uint64 i = 0; i < (1ULL << (1 + logN + lb + 2)); i ++)
        mpz_clear(beta[i]);
    for(uint64 i = 0; i < 4; i ++)
        mpz_clear(gamma[i]);
    free(poly_b);
    free(poly_d1);
    free(poly_d2);
    free(poly_r);
    free(R);
    free(Z);
    free(alpha1);
    free(alpha2);
    free(alpha3);
    free(tau);
    free(beta);
    free(gamma);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_dummy += time_diff;
	printf("%f ms dummy free \n", time_diff);
	timer = clock();

    printf("P(except commit part): %f ms\n", t_P);
    printf("V(except commit part): %f ms\n", t_V);
    return stat;
}


//Input: C1x(t) C1y(t); z
//Params: z, c
//Output: cC1x(t) cC1y(t); z
int sum_check_cipher_cmult(mpz_t rop, mpz_t *V, const mpz_t c, const mpz_t ri, const mpz_t z)
{
    int stat = 1;
    mpz_t inv, tmp1, pf;
    mpz_inits(inv, tmp1, pf, NULL);
    clock_t timer = clock();
    double time_diff, t_P = 0, t_V = 0;

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
    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover \n", time_diff);
	timer = clock();


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
            printf("Fail :: input not match %s %s\n", mpz_get_str(NULL, 16, ri), mpz_get_str(NULL ,16, tmp1));
        }
    }

    //return
    if(stat) {
        mpz_set(rop, pf);
    } else {
        mpz_set_ui(rop, 0);
    }
    mpz_clears(inv, tmp1, pf, NULL);

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier \n", time_diff);
	timer = clock();

    printf("P(except commit part): %f ms\n", t_P);
    printf("V(except commit part): %f ms\n", t_V);

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
    clock_t timer = clock();
    double time_diff, t_V = 0, t_P = 0;

    //Prover
    //update V :: from lsb
    mod_1neg(tmp1, z);
    for(int i = 0; i < 2; i ++){
        mod_mult(V[2 * i], V[2 * i], tmp1);
        mod_mult(V[2 * i + 1], V[2 * i + 1], z);
        mod_add(poly[i], V[2 * i], V[2 * i + 1]);
    }

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_P += time_diff;
	printf("%f ms Prover \n", time_diff);


    //Verifier
	timer = clock();
    mod_add(tmp1, poly[0], poly[1]);
    if(mpz_cmp(tmp1, ri)) {
        stat = 0;
        printf("Fail :: 1st sum check %s %s\n", mpz_get_str(NULL, 16, tmp1), mpz_get_str(NULL, 16, ri));
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

    time_diff = (double)(clock() - timer) / CLOCKS_PER_SEC * 1000;
    t_V += time_diff;
	printf("%f ms Verifier \n", time_diff);
	timer = clock();

    printf("P(except commit part): %f ms\n", t_P);
    printf("V(except commit part): %f ms\n", t_V);
   
    //free
    mpz_clears(tmp1, poly[0], poly[1], NULL);

    return stat;
}

//Input: Cx(t) Cy(t); z
//Params: logN, l1, l2, r[1 + logN + logQ], z[1 + logN + logQ]
//Output: [Cx](t) [Cy](t); z
int sum_check_cipher_rescale(mpz_t rop, mpz_t *C, const mpz_t ri, const int logN, const int l1, const int l2, const mpz_t *rn, const mpz_t *rd, const mpz_t *rb, const mpz_t val, const mpz_t *zn, const mpz_t *zdb,
	const mpz_t* commits, const uint64 lck, const uint64 key_num, const uint64 test_comm_num, const mpz_t* sks, const mpz_t gen)
{
    //ln = 1
    //ld = logN
    //lb = logQ
    return sum_check_poly_rounding(rop, C, ri, 1, logN, l1, l2, rn, rd, rb, val, zn, zdb, commits, lck, key_num, test_comm_num, sks, gen);
}


        for (uint64 j = 0; j < ((num * N * ubits) >> round); j ++) {
            mlmap_evaluation_N(c, 4, (num * N * ubits * 2) >> round, j, V_c);
            mlmap_evaluation_N(bb, 4, (num * N * ubits * 2) >> round, j, beta_b);
        }
        if (round < 1 + log_num) {
            for (uint64 j = 0; j < (num >> round); j ++) {
                mlmap_evaluation_N(bd, 3, (num * 2) >> round, j, beta_d);
                mlmap_evaluation_N(br, 3, (num * 2) >> round, j, beta_r);
            }
        } else if (round < 1 + log_num + logN) {
            for (uint64 j = 0; j < ((N / 2) >> (round - 1 - log_num)); j ++)
                mlmap_evaluation_N(t, 3, N >> (round - 1 - log_num), j, tau);
        } else {
            for (uint64 j = 0; j < ((ubits / 2) >> (round - 1 - log_num - logN)); j ++) {
                mlmap_evaluation_N(alpha1, 3, ubits >> (round - 1 - log_num - logN), j, a1);
                mlmap_evaluation_N(alpha2, 3, ubits >> (round - 1 - log_num - logN), j, a2);
            }
        }
        
*/
