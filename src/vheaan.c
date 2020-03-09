#include "vheaan.h"

int sum_check_verification(mpz_t rop, mpz_t **poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);

    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check\n%s\n%s\n%s\n", str, mpz_get_str(0, digit_rep, tmp), mpz_get_str(0, digit_rep, ri), mpz_get_str(0, digit_rep, PRIME));
        stat = 0;
    }
    
    if(stat) {
        polynomial_extrapolate_N(extrap_val, r[num_rounds - 1], poly[0], degree);
        for(int round = 1; round < num_rounds; round ++) {
            mod_add(tmp, poly[round][0], poly[round][1]);
            if(mpz_cmp(tmp, extrap_val)) {
                printf("Fail sumcheck %s :: %dth round %s %s\n", str, round, mpz_get_str(0, digit_rep, tmp), mpz_get_str(0, digit_rep, extrap_val));
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

int gkr_cipher_mult
(
    const mpz_t *V_0r, const mpz_t *V_0s, const mpz_t *V_0d,
    mpz_t *V_0c0, mpz_t *V_0c1,
    const mpz_t *V_1l, const mpz_t *V_1r, const mpz_t *V_1s, const mpz_t *V_1d,
    mpz_t *V_1c0, mpz_t * V_1c1,
    mpz_t *V_2, 

    const mpz_t *z_num, const mpz_t *z_N, const mpz_t *z_bits,
    const mpz_t *r_num, const mpz_t *r_N, const mpz_t *r_bits,
    const mpz_t P, const mpz_t EVK0, const mpz_t EVK1, 
    const int log_num, const int bits,
    const mpz_t val
)
{
    int num = 1 << log_num;
    digit_rep = 16;
    int log_bits = 0;
    while ((bits - 1) >> ++ log_bits);
    mpz_t *r = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 2 + logN + 1 + log_bits + 2));
    for (int i = 0; i < (int) (log_num + 1 + logN + 1 + log_bits + 2); i ++)
        mpz_init(r[i]);
    mpz_t v0rz, v0sz, v0dz, v1rz, v1sz, v1dz, v1lz, v2r, tmp1, tmp2, c00r, c01r, c10r, c11r, tmp0, tmp3, tmp4, tmp5, tmp6, BITMAX, BITINV;
    mpz_inits(v0rz, v0sz, v0dz, v1rz, v1sz, v1dz, v1lz, v2r, tmp1, tmp2, c00r, c01r, c10r, c11r, tmp0, tmp3, tmp4, tmp5, tmp6, BITMAX, BITINV, 0);
    mod_ui_pow_ui(BITMAX, 2, (1ULL << (log_bits + 2)) - 1);
    mod_inv_pow2(BITINV, bits);
    //First, get claimed value of V1*, V0* (w/o commits) on random point.
    //Pv-side (in the case of V_0r, Vf might should evaluate it too!)
    for (int i = 0; i < log_num; i ++)
        mpz_set(r[i], z_num[i + 1]);
    evaluate_V(v1rz, V_1r, log_num, r);
    evaluate_V(v1sz, V_1s, log_num, r);
    evaluate_V(v1dz, V_1d, log_num, r); 
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(r[i], z_num[i]);
    evaluate_V(v0sz, V_0s, log_num + 1, r);
    evaluate_V(v0dz, V_0d, log_num + 1, r);
    evaluate_V(v1lz, V_1l, log_num + 1, r);  
    evaluate_V(v0rz, V_0r, log_num + 1, r); //Vf also should do this.

    //Vf-side
    for (int i = 0; i < log_num + 2; i ++)
        mpz_set(r[i], r_num[i]);
    evaluate_V(v2r, V_2, log_num + 2, r);  

    //Commit Phase: this actually should be done at Fourth.
    //c00
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(r[i + logN + 1 + log_bits + 2], r_num[i + 1]);
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(r[i + log_bits + 2], r_N[i]);
    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(r[i], r_bits[i]);
    evaluate_V(c00r, V_0c0, log_num + 1 + logN + 1 + log_bits + 2, r);
    //c01
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(r[i + logN + log_bits + 2], r_num[i + 1]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(r[i + log_bits + 2], r_N[i]);
    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(r[i], r_bits[i]);
    evaluate_V(c01r, V_0c1, log_num + 1 + logN + log_bits + 2, r);
    //c10
    for (int i = 0; i < log_num; i ++)
        mpz_set(r[i + logN + 1 + log_bits + 2], r_num[i + 1]);
    for (int i = 0; i < (int) logN + 1; i ++)
        mpz_set(r[i + log_bits + 2], r_N[i]);
    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(r[i], r_bits[i]);
    evaluate_V(c10r, V_1c0, log_num + logN + 1 + log_bits + 2, r);
    //c11
    for (int i = 0; i < log_num; i ++)
        mpz_set(r[i + logN + log_bits + 2], r_num[i + 1]);
    for (int i = 0; i < (int) logN; i ++)
        mpz_set(r[i + log_bits + 2], r_N[i]);
    for (int i = 0; i < log_bits + 2; i ++)
        mpz_set(r[i], r_bits[i]);
    evaluate_V(c11r, V_1c1, log_num + logN + log_bits + 2, r);



    //Second, check the relationship between V1l, V1r and V0d.
    mod_1neg(tmp1, z_num[0]);
    mod_mult(tmp1, tmp1, EVK0);
    mod_mult(tmp2, z_num[0], EVK1); 
    mod_add(tmp1, tmp1, tmp2);
    mod_mult(tmp1, tmp1, v1rz); //tmp1 <- ((1-r)*EVK0+r*EVK1)*v1rr
    mod_mult(tmp2, P, v1lz); //tmp2 <- P*v1lr
    mod_add(tmp1, tmp1, tmp2);
    if (mpz_cmp(tmp1, v0dz)) {
        printf("F.A.I.L.\n%s\n%s\n", mpz_get_str(0, digit_rep, tmp1), mpz_get_str(0, digit_rep, v0dz));
    } else {
        printf("S.U.C.C.E.E..ee.s.:-)\n");
    }
   
    
    //Third, check the GKR protocols connecting lower layers to upper layers.
    // Check 14 GKR protocols simultaneously.
    mpz_t *beta_num = (mpz_t *) malloc(sizeof(mpz_t) * num * 2),
          *beta_num_2 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *beta_N = (mpz_t *) malloc(sizeof(mpz_t) * N * 2),
          *beta_bits = (mpz_t *) malloc(sizeof(mpz_t) * (1ULL << log_bits) * 4),
          *alpha_4q = (mpz_t *) malloc(sizeof(mpz_t) * (1ULL << log_bits) * 4),
          *alpha_2q = (mpz_t *) malloc(sizeof(mpz_t) * (1ULL << log_bits) * 4),
          *alpha_q = (mpz_t *) malloc(sizeof(mpz_t) * (1ULL << log_bits) * 4),
          *tau_N = (mpz_t *) malloc(sizeof(mpz_t) * N),
          *tau_2N = (mpz_t *) malloc(sizeof(mpz_t) * N * 2),
          *V_2_00 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_2_01 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_2_10 = (mpz_t *) malloc(sizeof(mpz_t) * num),
          *V_2_11 = (mpz_t *) malloc(sizeof(mpz_t) * num);
    mpz_t **poly_1l = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 2)),
          **poly_1d2 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 2)),
          **poly_1d = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
          **poly_1s0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
          **poly_1b0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
          **poly_1r = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
          **poly_1s1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
          **poly_1b1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
          **poly_0d = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
          **poly_0s0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
          **poly_0b0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
          **poly_0r = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2)),
          **poly_0s1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2)),
          **poly_0b1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2));
    mpz_t bn[4], bn2[4], bN[4], bb[4], a4q[4], a2q[4], aq[4], tN1[4], t2N[4], tN2[4], c00[4], c01[4], c10[4], c11[4], v200[4], v201[4], v210[4], v211[4];
printf("1\n");
    for (int i = 0; i < num; i ++) {
        mpz_inits(V_2_00[i], V_2_01[i], V_2_10[i], V_2_11[i], beta_num[i], 0);
        mpz_set(V_2_00[i], V_2[i * 4]);
        mpz_set(V_2_01[i], V_2[i * 4 + 1]);
        mpz_set(V_2_10[i], V_2[i * 4 + 2]);
        mpz_set(V_2_11[i], V_2[i * 4 + 3]);
    }
    for (int i = 0; i < num * 2; i ++)
        mpz_init(beta_num[i]);
    for (int i = 0; i < num; i ++)
        mpz_init(beta_num_2[i]);
    for (int i = 0; i < (int) N; i ++)
        mpz_inits(tau_N[i], 0);
    for (int i = 0; i < (int) N * 2; i ++)
        mpz_inits(beta_N[i], tau_2N[i], 0);
    for (int i = 0; i < (1 << (log_bits + 2)); i ++)
        mpz_inits(beta_bits[i], alpha_4q[i], alpha_2q[i], alpha_q[i], 0);
    initialize_beta(beta_num, log_num + 1, z_num);
    for (int i = 0; i < log_num; i ++)
        mpz_set(r[i], z_num[i + 1]);
    initialize_beta(beta_num_2, log_num, r);
    initialize_beta(beta_N, logN + 1, z_N);
    initialize_tau(tau_N, logN, val);
    initialize_tau(tau_2N, logN + 1, val);
    initialize_beta(beta_bits, log_bits + 2, z_bits);
    initialize_alpha(alpha_4q, log_bits + 2, bits * 4);
    initialize_alpha(alpha_2q, log_bits + 2, bits * 2);
    initialize_alpha(alpha_q, log_bits + 2, bits);
    for (int i = 0; i < 4; i ++)
        mpz_inits(bn[i], bn2[i], bN[i], bb[i], a4q[i], a2q[i], aq[i], tN1[i], t2N[i], tN2[i], c00[i], c01[i], c10[i], c11[i], v210[i], v200[i], v201[i], v211[i], 0);
    for (int i = 0; i < log_num + 2; i ++) {
        poly_1l[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_1d2[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_1l[i][j], 0);
            mpz_init_set_ui(poly_1d2[i][j], 0);
        }
    }
    for (int i = 0; i < (int) (log_num + logN + 1 + log_bits + 2); i ++) {
        poly_1d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_1s0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_1b0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_1d[i][j], 0);
            mpz_init_set_ui(poly_1s0[i][j], 0);
            mpz_init_set_ui(poly_1b0[i][j], 0);
        }
    }
    for (int i = 0; i < (int) (log_num + logN + log_bits + 2); i ++) {
        poly_1r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_1s1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_1b1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_1r[i][j], 0);
            mpz_init_set_ui(poly_1s1[i][j], 0);
            mpz_init_set_ui(poly_1b1[i][j], 0);
        }
    }
    for (int i = 0; i < (int) (log_num + 1 + logN + 1 + log_bits + 2); i ++) {
        poly_0d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_0s0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_0b0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_0d[i][j], 0);
            mpz_init_set_ui(poly_0s0[i][j], 0);
            mpz_init_set_ui(poly_0b0[i][j], 0);
        }
    }
    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits + 2); i ++) {
        poly_0r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_0s1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        poly_0b1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
        for (int j = 0; j < 4; j ++) {
            mpz_init_set_ui(poly_0r[i][j], 0);
            mpz_init_set_ui(poly_0s1[i][j], 0);
            mpz_init_set_ui(poly_0b1[i][j], 0);
            mpz_inits(poly_0r[i][j], poly_0s1[i][j], poly_0b1[i][j], 0);
        }
    }
    printf("memory allocated\n");

    uint64 num_terms_num = 1ULL << (log_num + 1);
    uint64 num_terms_N = 1ULL << (logN + 1);
    uint64 num_terms_bits = 1ULL << (log_bits + 2);
    for (int round = 0; round < log_num + 1; round ++) {
        for (uint64 j = 0; j < num_terms_num / 2; j ++) {
            //evaluate mle
            mlmap_evaluation_N(bn, 4, num_terms_num, j, beta_num);
            if ((round < log_num) && (!(j & 1)))
                mlmap_evaluation_N(bn2, 4, num_terms_num >> 1, j >> 1, beta_num_2);

            for (uint64 k = 0; k < num_terms_N; k ++) {
                //evaluate mle
                mlmap_evaluation_N(bN, 4, 0, k, beta_N);
                mlmap_evaluation_N(tN1, 4, 0, k >> 1, tau_N);
                mlmap_evaluation_N(tN2, 4, 0, k & (num_terms_N / 2 - 1), tau_2N);
                mlmap_evaluation_N(t2N, 4, 0, k, tau_2N);


                for (uint64 t = 0; t < num_terms_bits; t ++) {
                    //evaluate mle
                    mlmap_evaluation_N(a4q, 4, 0, t, alpha_4q);
                    mlmap_evaluation_N(a2q, 4, 0, t, alpha_2q);
                    mlmap_evaluation_N(aq, 4, 0, t, alpha_q);
                    mlmap_evaluation_N(bb, 4, 0, t, beta_bits);
                    
                    if ((round < log_num) && (!(j & 1))) {//앞쪽 randomness를 활용
                        mlmap_evaluation_N(c10, 4, num_terms_num << (-1 + logN + 1 + log_bits + 2), ((((j >> 1) << (logN + 1)) + k) << (log_bits + 2)) + t, V_1c0);

                        
                        for (uint64 l = 0; l < 4; l ++) {
                            //poly_1b0
                            mod_sub_ui(tmp1, c10[l], 1);
                            mod_mult(tmp1, tmp1, c10[l]);
                            mod_mult(tmp1, tmp1, bn2[l]);
                            mod_mult(tmp1, tmp1, bN[l]);
                            mod_mult(tmp1, tmp1, bb[l]);
                            mod_add(poly_1b0[round][l], poly_1b0[round][l], tmp1);

                            //poly_1d
                            mod_mult(tmp1, bn2[l], t2N[l]);
                            mod_mult(tmp1, tmp1, a4q[l]);
                            mod_mult(tmp1, tmp1, c10[l]);
                            mod_add(poly_1d[round][l], poly_1d[round][l], tmp1);

                            //poly_1s0
                            if (k >> logN) {
                                mod_ui_sub(tmp1, 0, c00[l]);
                                mod_mult(tmp1, tmp1, a4q[l]);
                                mod_mult(tmp1, tmp1, tN2[l]);
                                mod_mult(tmp1, tmp1, bn2[l]);
                                if (!t){
                                    mod_mult(tmp2, bn2[l], tN2[l]);
                                    mod_mult(tmp2, tmp2, BITMAX);
                                    mod_add(tmp1, tmp1, tmp2);
                                }
                            } else {
                                mod_mult(tmp1, a4q[l], tN2[l]);
                                mod_mult(tmp1, tmp1, c00[l]);
                                mod_mult(tmp1, tmp1, bn2[l]);
                            }
                            mod_add(poly_1s0[round][l], poly_1s0[round][l], tmp1);
                        }

                        if (!(k & 1)) { //얘도 그냥 앞쪽을 활용하자: 이를 통해 bN을 쉽게 계산..
                            mlmap_evaluation_N(c11, 4, num_terms_num << (-1 + logN + log_bits + 2), ((((j >> 1) << logN) + (k >> 1)) << (log_bits + 2)) + t, V_1c1);
                            for (uint64 l = 0; l < 4; l ++) {
                                //poly_1b1
                                mod_sub_ui(tmp1, c11[l], 1);
                                mod_mult(tmp1, tmp1, c11[l]);
                                mod_mult(tmp1, tmp1, bn2[l]);
                                mod_mult(tmp1, tmp1, bN[l]);
                                mod_mult(tmp1, tmp1, bb[l]);
                                mod_add(poly_1b1[round][l], poly_1b1[round][l], tmp1);

                                //poly_1s1
                                mod_mult(tmp1, a4q[l], bn2[l]);
                                mod_mult(tmp1, tmp1, tN1[l]);
                                mod_mult(tmp1, tmp1, c11[l]);
                                mod_add(poly_1s1[round][l], poly_1s1[round][l], tmp1);

                                //poly_1r
                                mod_sub(tmp1, aq[l], bn2[l]);
                                mod_mult(tmp1, tmp1, tN1[l]);
                                mod_mult(tmp1, tmp1, c11[l]);
                                mod_add(poly_1r[round][l], poly_1r[round][l], tmp1);
                            }

                        }
                    }


                    mlmap_evaluation_N(c00, 4, num_terms_num << (logN + 1 + log_bits + 2), (((j << (logN + 1)) + k) << (log_bits + 2)) + t, V_0c0);

                    
                    for (uint64 l = 0; l < 4; l ++) {
                        //poly_0b0
                        mod_sub_ui(tmp1, c00[l], 1);
                        mod_mult(tmp1, tmp1, c00[l]);
                        mod_mult(tmp1, tmp1, bn[l]);
                        mod_mult(tmp1, tmp1, bN[l]);
                        mod_mult(tmp1, tmp1, bb[l]);
                        mod_add(poly_0b0[round][l], poly_0b0[round][l], tmp1);

                        //poly_0d
                        mod_mult(tmp1, bn[l], t2N[l]);
                        mod_mult(tmp1, tmp1, a4q[l]);
                        mod_mult(tmp1, tmp1, c00[l]);
                        mod_add(poly_0d[round][l], poly_0d[round][l], tmp1);

                        //poly_0s0
                        if (k >> logN) {
                            mod_sub(tmp1, PRIME, bn[l]);
                            mod_mult(tmp1, tmp1, tN2[l]);
                            mod_mult(tmp1, tmp1, a4q[l]);
                            mod_mult(tmp1, tmp1, c00[l]);
                            if (!t) {
                                mod_mult(tmp2, bn[l], tN2[l]);
                                mod_mult(tmp2, tmp2, BITMAX);
                                mod_add(tmp1, tmp1, tmp2);
                            }
                        } else {
                            mod_mult(tmp1, a4q[l], tN2[l]);
                            mod_mult(tmp1, tmp1, c00[l]);
                            mod_mult(tmp1, tmp1, bn[l]);
                        }
                        mod_add(poly_0s0[round][l], poly_0s0[round][l], tmp1);
                    }

                    if (!(k & 1)) { //얘도 그냥 앞쪽을 활용하자: 이를 통해 bN을 쉽게 계산..
                        mlmap_evaluation_N(c01, 4, num_terms_num << (logN + log_bits + 2), (((j << logN) + (k >> 1)) << (log_bits + 2)) + t, V_0c1);
                        for (uint64 l = 0; l < 4; l ++) {
                            //poly_0b1
                            mod_sub_ui(tmp1, c01[l], 1);
                            mod_mult(tmp1, tmp1, c01[l]);
                            mod_mult(tmp1, tmp1, bn[l]);
                            mod_mult(tmp1, tmp1, bN[l]);
                            mod_mult(tmp1, tmp1, bb[l]);
                            mod_add(poly_0b1[round][l], poly_0b1[round][l], tmp1);

                            //poly_0s1
                            mod_mult(tmp1, a4q[l], bn[l]);
                            mod_mult(tmp1, tmp1, tN1[l]);
                            mod_mult(tmp1, tmp1, c01[l]);
                            mod_add(poly_0s1[round][l], poly_0s1[round][l], tmp1);

                            //poly_0r
                            mod_sub(tmp1, a2q[l], aq[l]);
                            mod_mult(tmp1, tmp1, BITINV);
                            mod_mult(tmp1, tmp1, bn[l]);
                            mod_mult(tmp1, tmp1, tN1[l]);
                            mod_mult(tmp1, tmp1, c01[l]);
                            mod_add(poly_0r[round][l], poly_0r[round][l], tmp1);
                        }
                    }
                }
            }

            if ((round < log_num) && (j & 1)) {
                mlmap_evaluation_N(v200, 4, num_terms_num >> 1, j >> 1, V_2_00);
                mlmap_evaluation_N(v201, 4, num_terms_num >> 1, j >> 1, V_2_01);
                mlmap_evaluation_N(v210, 4, num_terms_num >> 1, j >> 1, V_2_10);
                mlmap_evaluation_N(v211, 4, num_terms_num >> 1, j >> 1, V_2_11);

                for (uint64 l = 0; l < 4; l ++) {
                    //poly_1d
                    mod_mult(tmp1, v201[l], v211[l]);
                    mod_mult(tmp1, tmp1, bn[l]);
                    mod_add(poly_1d[round][l], poly_1d[round][l], tmp1);

                    //poly_1l
                    mod_1neg(tmp1, z_num[log_num]);
                    mod_mult(tmp1, tmp1, v200[l]);
                    mod_mult(tmp1, tmp1, v210[l]);
                    mod_mult(tmp2, v211[l], v200[l]);
                    mod_mult(tmp3, v201[l], v210[l]);
                    mod_add(tmp2, tmp3, tmp2);
                    mod_mult(tmp2, tmp2, z_num[log_num]);
                    mod_add(tmp1, tmp1, tmp2);
                    mod_add(poly_1l[round][l], poly_1l[round][l], tmp1);
                }
            }
        }

        num_terms_num >>= 1;
        update_V(beta_num, num_terms_num, r_num[log_num + 2 - round - 1]);

        if (round < log_num) {
            update_V(beta_num_2, num_terms_num >> 1, r_num[log_num + 2 - round - 1]);

            update_V(V_2_00, num_terms_num >> 1, r_num[log_num + 2 - round - 1]);
            update_V(V_2_01, num_terms_num >> 1, r_num[log_num + 2 - round - 1]);
            update_V(V_2_10, num_terms_num >> 1, r_num[log_num + 2 - round - 1]);
            update_V(V_2_11, num_terms_num >> 1, r_num[log_num + 2 - round - 1]);

            update_V(V_1c0, num_terms_num << (-1 + logN + 1 + log_bits + 2), r_num[log_num + 2 - round - 1]);
            update_V(V_1c1, num_terms_num << (-1 + logN + log_bits + 2), r_num[log_num + 2 - round - 1]);
        }
        update_V(V_0c0, num_terms_num << (logN + 1 + log_bits + 2), r_num[log_num + 2 - round - 1]);
        update_V(V_0c1, num_terms_num << (logN + log_bits + 2), r_num[log_num + 2 - round - 1]);
    }
    //Verify
    int stat = 1;
    mpz_set_ui(tmp1, 0);
    for (int i = 0; i < log_num + 1; i ++)
        mpz_set(r[i], r_num[i + 1]);
    stat = stat ? sum_check_verification(tmp2, poly_0b0, tmp1, 4, log_num + 1, r, "V0b0") : 0;
    if(stat)
        printf("V0b0 Verified\n");
    stat = stat ? sum_check_verification(tmp2, poly_0b1, tmp1, 4, log_num + 1, r, "V0b1") : 0;
    if(stat)
        printf("V0b1 Verified\n");
    mpz_set(tmp1, v0dz);
    stat = stat ? sum_check_verification(tmp2, poly_0d, tmp1, 4, log_num + 1, r, "V0d") : 0;
    if(stat)
        printf("V0d Verified\n");
    mpz_set(tmp1, v0rz);
    stat = stat ? sum_check_verification(tmp2, poly_0r, tmp1, 4, log_num + 1, r, "V0r") : 0;
    if(stat)
        printf("V0r Verified\n");
    stat = 1;
    mpz_set(tmp1, v0sz);
    stat = stat ? sum_check_verification(tmp2, poly_0s0, tmp1, 4, log_num + 1, r, "V0s0") : 0;
    if(stat)
        printf("V0s0 Verified\n");
    stat = 1;
    stat = stat ? sum_check_verification(tmp2, poly_0s1, tmp1, 4, log_num + 1, r, "V0s1") : 0;
    if(stat)
        printf("V0s1 Verified\n");

    stat = 1;
    mpz_set_ui(tmp1, 0); 
    for (int i = 0; i < log_num; i ++)
        mpz_set(r[i], r_num[i + 2]);
    stat = stat ? sum_check_verification(tmp2, poly_1b0, tmp1, 4, log_num, r, "V1b0") : 0;
    if(stat)
        printf("V1b0 Verified\n");
    stat = 1;
    stat = stat ? sum_check_verification(tmp2, poly_1b1, tmp1, 4, log_num, r, "V1b1") : 0;
    if(stat)
        printf("V1b1 Verified\n");
    stat = 1;
    mpz_set(tmp1, v0dz);
    stat = stat ? sum_check_verification(tmp2, poly_1d, tmp1, 4, log_num, r, "V1d") : 0;
    if(stat)
        printf("V1d Verified\n");
    mpz_set(tmp1, v0rz);
    stat = stat ? sum_check_verification(tmp2, poly_1r, tmp1, 4, log_num, r, "V1r") : 0;
    if(stat)
        printf("V1r Verified\n");
    stat = 1;
    mpz_set(tmp1, v0sz);
    stat = stat ? sum_check_verification(tmp2, poly_1s0, tmp1, 4, log_num, r, "V1s0") : 0;
    if(stat)
        printf("V1s0 Verified\n");
    stat = 1;
    stat = stat ? sum_check_verification(tmp2, poly_1s1, tmp1, 4, log_num, r, "V1s1") : 0;
    if(stat)
        printf("V1s1 Verified\n");




    for (int round = 0; round < (int) logN + 1; round ++) {
        for (uint64 j = 0; j < num_terms_N / 2; j ++) {
            for (uint64 k = 0; k < num_terms_bits / 2; k ++) {
                
            }
        }        
    }
    for (int round = 0; round < log_bits + 2; round ++) {
        for (uint64 j = 0; j < num_terms_bits / 2; j ++) {

        }
    }

    
    //Fourth, check the clamied value of the commited MLE is consistent with the committed values.
    //Currently, ommitted.

    
    //Fifth, Check the result of GKR protocol w/ the precalculated value.
    if (mpz_cmp(tmp5, v2r)) {
        printf("F.A.I.L. (V2)\n%s\n%s\n", mpz_get_str(0, digit_rep, tmp5), mpz_get_str(0, digit_rep, v2r));
    } else {
        printf("Succececesesese\n");
    }

    return 0;
}

//    //Fifth, get the claimed value of V2[0]~V2[3].
//    //!can boost up here!
//    for (int i = 0; i < log_num; i ++)
//        mpz_set(r[i + 2], r_num[i + 1]);
//    for (int i = 0; i < 2; i ++)
//        mpz_set_ui(r[i], (0 >> i) & 1);
//    evaluate_V(tmp0, V_2, log_num + 2, r);
//    for (int i = 0; i < 2; i ++)
//        mpz_set_ui(r[i], (1 >> i) & 1);
//    evaluate_V(tmp1, V_2, log_num + 2, r);
//    for (int i = 0; i < 2; i ++)
//        mpz_set_ui(r[i], (2 >> i) & 1);
//    evaluate_V(tmp2, V_2, log_num + 2, r);
//    for (int i = 0; i < 2; i ++)
//        mpz_set_ui(r[i], (3 >> i) & 1);
//    evaluate_V(tmp3, V_2, log_num + 2, r);
//
//
//    //Check the consistency between the claimed V2[i] and V1d.
//    mpz_t *poly_b = (mpz_t *) malloc(sizeof(mpz_t) * ());
//    mod_mult(tmp4, tmp1, tmp3);
//    printf("Hello~~~~~\n%s\n%s\n",
//        mpz_get_str(0, digit_rep, tmp4),
//        mpz_get_str(0, digit_rep, v1dr));
//    mod_mult(tmp5, V_2[1], V_2[3]);
//    printf("YAAAAA\n%s\n%s\n", 
//        mpz_get_str(0, digit_rep, V_1d[0]),
//        mpz_get_str(0, digit_rep, tmp5));
//
//
//    //Check the consistency between the claimed V2[i] and V1l.
//    mod_1neg(tmp4, r_num[0]);
//    mod_mult(tmp4, tmp4, tmp0);
//    mod_mult(tmp4, tmp4, tmp2);//(1-r)v2(0)v2(2)
//    mod_mult(tmp5, tmp1, tmp2);
//    mod_mult(tmp6, tmp0, tmp3);
//    mod_add(tmp5, tmp5, tmp6);
//    mod_mult(tmp5, tmp5, r_num[0]);//r(v2(0)v2(3)+v2(1)v2(2))
//    mod_add(tmp4, tmp4, tmp5);
//    if (mpz_cmp(tmp4, v1lr)) {
//        printf("F.A.I.L. (v1l)\n%s\n%s\n", mpz_get_str(0, digit_rep, tmp4), mpz_get_str(0, digit_rep, v1lr));
//    } else {
//        printf("S.u.cc..c.e.s.ssss\n");
//    }
//
//
//
//    //Evaluate V2[r] from V2[i]
//    mod_1neg(tmp4, r_2[1]);
//    mod_1neg(tmp5, r_2[0]);
//    mod_mult(tmp4, tmp4, tmp5);
//    mod_mult(tmp5, tmp4, tmp0);//(1-r1)(1-r0)V2(00)
//    mod_1neg(tmp4, r_2[1]);
//    mod_mult(tmp4, tmp4, r_2[0]);
//    mod_mult(tmp4, tmp4, tmp1);//(1-r1)r0V2(01)
//    mod_add(tmp5, tmp5, tmp4);
//    mod_1neg(tmp4, r_2[0]);
//    mod_mult(tmp4, tmp4, r_2[1]);
//    mod_mult(tmp4, tmp4, tmp2);//r1(1-r0)V2(10)
//    mod_add(tmp5, tmp5, tmp4);
//    mod_mult(tmp4, r_2[0], r_2[1]);
//    mod_mult(tmp4, tmp4, tmp3);//r1r0V2(11)
//    mod_add(tmp5, tmp5, tmp4);
 
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
          *tmps2 = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num + logN + log_bits));
          //*r_d = (mpz_t *) malloc(sizeof(mpz_t) * (1 + log_num));

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
    printf("GKR Rescale: Proof generated\n");

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
    printf("GKR Rescale: Verified (%d)\n", stat);

    return stat;
}
//
//int gkr_cipher_mult
//(
//    const mpz_t *V_0r, const mpz_t *V_0s, const mpz_t *V_0d,
//    const mpz_t *V_0c0, const mpz_t *V_0c1,
//    const mpz_t *V_1l, const mpz_t *V_1r, const mpz_t *V_1s, const mpz_t *V_1d,
//    const mpz_t *V_1c0, const mpz_t * V_1c1,
//    const mpz_t *V_2, 
//
//    const mpz_t *r_num, 
//    const mpz_t *r_0, const mpz_t *r_1, const mpz_t *r_2,
//    const mpz_t P, const mpz_t EVK1, const mpz_t EVK2, 
//    const int log_num, const int bits,
//    const mpz_t val
//)
//{
//    digit_rep = 16;
//    int log_bits = 0;
//    while ((bits - 1) >> ++ log_bits);
//    int ubits = 1 << log_bits, num = 1 << log_num, stat = 1;
//    uint64 num_terms;
//
//    mpz_t tmp1, tmp2, scale, c00r, c01r, c10r, c11r, a1r, a2r, bhr, btr, tr, v2r, v0rr;
//    mpz_inits(tmp1, tmp2, scale, c00r, c01r, c10r, c11r, a1r, a2r, bhr, btr, tr, v2r, v0rr, 0);
//
//    mpz_t bh[4], bt[4], c00[4], c01[4], c10[4], c11[4], a1[3], a2[3], t[3];
//    for (int i = 0; i < 4; i ++)
//        mpz_inits(bb[i], bt[i], c00[i], c01[i], c10[i], c11[i], 0);
//    for (int i = 0; i < 3; i ++)
//        mpz_inits(bh[i], a1[i], a2[i], t[i], 0);
//    
//    mpz_t *beta_Al = (mpz_t *) malloc(sizeof(mpz_t) * (num *2)),
//          *beta_As = (mpz_t *) malloc(sizeof(mpz_t) * num),
//          *beta_Bl = (mpz_t *) malloc(sizeof(mpz_t) * (N * 2)),
//          *beta_Bs = (mpz_t *) malloc(sizeof(mpz_t) * N),
//          *beta_C = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
//
//          *alpha1 = (mpz_t *) malloc(sizeof(mpz_t) * ubtis * 4),
//          *alpha2 = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
//          *alpha3 = (mpz_t *) malloc(sizeof(mpz_t) * ubits * 4),
//
//          *tau_l0 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
//          *tau_l1 = (mpz_t *) malloc(sizeof(mpz_t) * 2 * N),
//          *tau_s = (mpz_t *) malloc(sizeof(mpz_t) * N),
//
//          **poly_1d = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
//          **poly_1s0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
//          **poly_1s1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
//          **poly_1r = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
//          **poly_0d = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
//          **poly_0s0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
//          **poly_0s1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2)),
//          **poly_0r = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2)),
//          **poly_1b0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + 1 + log_bits + 2)),
//          **poly_1b1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + logN + log_bits + 2)),
//          **poly_0b0 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + 1 + log_bits + 2)),
//          **poly_0b1 = (mpz_t **) malloc(sizeof(mpz_t *) * (log_num + 1 + logN + log_bits + 2)),
//          
//          *tmps1 = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + 1 + log_bits + 2)),
//          *tmps2 = (mpz_t *) malloc(sizeof(mpz_t) * (log_num + 1 + logN + 1 + log_bits + 2));
//
//    for (int i = 0; i < 2 * num; i ++) 
//        mpz_init(beta_Al[i]);
//    for (int i = 0; i < num; i ++)
//        mpz_init(beta_As[i]);
//    for (int i = 0; i < N * 2; i ++)
//        mpz_init(beta_Bl[i]);
//    for (int i = 0; i < N; i ++)
//        mpz_init(beta_Bs[i]);
//    for (int i = 0; i < ubits * 4; i ++)
//        mpz_init(beta_C[i]);
//    for (int i = 0; i < ubits * 4; i ++)
//        mpz_inits(alpha1[i], alpha2[i], alpha3[i], 0);
//    for (int i = 0; i < 2 * N; i ++)
//        mpz_inits(tau_l0[i], tau_l1[i], 0);
//    for (int i = 0; i < N; i ++)
//        mpz_init(tau_s[i]);
//    for (int i = 0; i < (int) (log_num + logN + 1 + log_bits + 2); i ++) {
//        poly_1d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_1s0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_1b0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
//        for (int j = 0; j < 3; j ++)
//            mpz_inits(poly_1d[i][j], poly_1s0[i], 0);
//        for (int j = 0; j < 4; j ++)
//            mpz_init(poly_1b0[i][j]);
//    }
//    for (int i = 0; i < (int) (log_num + logN + log_bits + 2); i ++) {
//        poly_1r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_1s1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_1b1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
//        for (int j = 0; j < 3; j ++)
//            mpz_inits(poly_1r[i][j], poly_1s1[i], 0);
//        for (int j = 0; j < 4; j ++)
//            mpz_init(poly_1b1[i][j]);
//    }
//    for (int i = 0; i < (int) (log_num + 1 + logN + 1 + log_bits + 2); i ++) {
//        poly_0d[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_0s0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_0b0[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
//        for (int j = 0; j < 3; j ++)
//            mpz_inits(poly_0d[i][j], poly_0s0[i], 0);
//        for (int j = 0; j < 4; j ++)
//            mpz_init(poly_0b0[i][j]);
//    }
//    for (int i = 0; i < (int) (log_num + 1 + logN + log_bits + 2); i ++) {
//        poly_0r[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_0s1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 3);
//        poly_0b1[i] = (mpz_t *) malloc(sizeof(mpz_t) * 4);
//        for (int j = 0; j < 3; j ++)
//            mpz_inits(poly_0r[i][j], poly_0s1[i], 0);
//        for (int j = 0; j < 4; j ++)
//            mpz_init(poly_0b1[i][j]);
//    }
//    for (int i = 0; i < (int) (log_num + 1 + logN + 1 + log_bits + 2); i ++)
//        mpz_inits(tmps1[i], tmps2[i], 0);
//    printf("GKR Mult: Memory allocated\n");
//
//
//    //init mle
//    initialize_beta(beta_Al, log_num + 1, );
//    initialize_beta(beta_As, log_num, );
//    initialize_beta(beta_Bl, logN + 1, );
//    initialize_beta(beta_Bs, logN, );
//    initialize_beta(beta_C, log_bits + 2, );
//    initialize_alpha(alpha1, log_bits + 2, bits);
//    initialize_alpha(alpha2, log_bits + 2, 2 * bits);
//    initialize_alpha(alpha3, log_bits + 2, 4 * bits);
//    initialize_tau(tau_l0, logN + 1, val)
//    initialize_tau2(tau_l1, logN + 1, val)
//    initialize_tau(tau_s, logN, val)
//    for (int i = 0; i < (int) N * num * 2 * (1 << log_bits); i ++)
//        mpz_init_set(V_c[i], V_c_in[i]);
//    mod_inv_pow2(scale, dif_bits);
//    printf("GKR Mult: MLE initialized\n");
//
//
//    //eval mle
////    evaluate_V(v2r, V_2, log_num + 2, );
////    evaluate_V(v0rr, V_0r, log_num + 1, );
////
////    evaluate_V(c00r, V_0c0, log_num + logN + 1 + log_bits + 2, );
////    evaluate_V(c01r, V_0c1, log_num + logN + log_bits + 2, );
////    evaluate_V(c10r, V_1c0, log_num + 1 + logN + 1 + log_bits + 2, );
////    evaluate_V(c11r, V_1c1, log_num + 1 + logN + log_bits + 2, );
////
////    evaluate_beta(bhr, r_b, r_c, log_num);
////    evaluate_beta(btr, r_b, r_c, logN + log_bits + 2);
////
////    evaluate_tau(tr, val, , logN);
////    
////    evaluate_alpha(a1r, bits, , log_bits + 2);
////    evaluate_alpha(a2r, 2 * bits, , log_bits + 2);
//    printf("GKR Mult: MLE Precomputed\n");
//
//
//
//    num_terms = 1ULL << (log_num + 1 + logN + 1 + log_bits + 2);
//    for (uint64 j = 0; j < num_terms / 2; j ++) {
//        //linear computation
//        mlmap_evaluation_N(bAl, 4, num_terms >> (logN + 1 + log_bits + 2), j >> (logN + 1 + log_bits + 2), beta_Al);
//        mlmap_evaluation_N(bAs, 4, 0, (j >> (logN + 1 + log_bits + 2)) & (num * 2 - 1), beta_Als);
//        mlmap_evaluation_N(bBl, 4, 0, (j >> (log_bits + 2)) & (N * 4 - 1), beta_Bl);
//        mlmap_evaluation_N(bBs, 4, 0, (j >> (log_bits + 2)) & (N * 2 - 1), beta_Bs);
//        mlmap_evaluation_N(bC, 4, 0, j & (8 * ubits - 1), beta_C);
//
//        mlmap_evaluation_N(c00, 4, num_terms, j, V_0c0);
//        mlmap_evaluation_N(c01, 4, num_terms, j, V_0c1);
//
//        mlmap_evaluation_N(a1, 3, 0, j & (8 * ubits - 1), alpha1);
//        mlmap_evaluation_N(a2, 3, 0, j & (8 * ubits - 1), alpha2);
//        mlmap_evaluation_N(a3, 3, 0, j & (8 * ubits - 1), alpha3);
//
//        mlmap_evaluation_N(tl0, 3, 0, (j >> (log_bits + 2)) & (N * 4 - 1), tau_l0);
//        mlmap_evaluation_N(tl1, 3, 0, (j >> (log_bits + 2)) & (N * 4 - 1), tau_l1);
//        mlmap_evaluation_N(ts, 3, 0, (j >> (log_bits + 2)) & (N * 2 - 1), tau_s);
//
//
//        //update poly
//        for (int k = 0; k < 3; k ++) {
//            mod_mult(tmp1, bAl[k], a3[k]);
//            mod_mult(tmp1, tmp1, bBl[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, tl0[k]);
//            mod_mult(tmp1, tmp1, c00[k]);
//            mod_add(poly_0d[k], poly_0d[k], tmp1);
//        }
//        for (int k = 0; k < 3; k ++) {
//            mod_mult(tmp1, bAl[k], a3[k]);
//            mod_mult(tmp1, tmp1, bBl[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, tl1[k]);
//            mod_mult(tmp1, tmp1, c00[k]);
//            mod_add(poly_0s0[k], poly_0s0[k], tmp1);
//        }
//        for (int k = 0; k < 3; k ++) {
//            mod_mult(tmp1, bAl[k], a3[k]);
//            mod_mult(tmp1, tmp1, bBs[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, ts[k]);
//            mod_mult(tmp1, tmp1, c01[k]);
//            mod_add(poly_0s1[k], poly_0s1[k], tmp1);
//        }
//        for (int k = 0; k < 3; k ++) {
//            mod_sub(tmp1, a2[k], a1[k]);
//            mod_add(tmp1, tmp1, P);
//            mod_mult(tmp1, tmp1, bAl[k]);
//            mod_mult(tmp1, tmp1, bBs[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, a3[k]);
//            mod_mult(tmp1, tmp1, ts[k]);
//            mod_mult(tmp1, tmp1, c01[k]);
//            mod_add(poly_0r[k], poly_0r[k], tmp1);
//        }
//        for (int k = 0; k < 4; k ++) {
//            mod_sub_ui(tmp1, c00[k], 1);
//            mod_mult(tmp1, tmp1, c00[k]);
//            mod_mult(tmp1, tmp1, bAl[k]);
//            mod_mult(tmp1, tmp1, bBl[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, tl[k]);
//            mod_add(poly_0b0[k], poly_0b0[k], tmp1);
//        }
//        for (int k = 0; k < 4; k ++) {
//            mod_sub_ui(tmp1, c01[k], 1);
//            mod_mult(tmp1, tmp1, c01[k]);
//            mod_mult(tmp1, tmp1, bAl[k]);
//            mod_mult(tmp1, tmp1, bBs[k]);
//            mod_mult(tmp1, tmp1, bC[k]);
//            mod_mult(tmp1, tmp1, ts[k]);
//            mod_add(poly_0b1[k], poly_0b1[k], tmp1);
//        }
//    }
//
//    //update mle
//    num_terms >>= 1;
//    mpz_set(r, r_1[log_num + 1 + logN + 1 + log_bits + 2 - 1]);
//    update_V(beta_Al, num_terms >> (logN + 1 + log_bits + 2), r);
//    update_V(V_1c0, num_terms, r);
//    update_V(V_1c1, num_terms, r);
//
//
//    for (int round = 0; round < log_num; i ++) {
//        for (uint64 j = 0; j < num_terms / 2; j ++) {
//            //linear computation
//            mlmap_evaluation_N(bAl, 4, num_terms >> (logN + 1 + log_bits + 2), j >> (logN + 1 + log_bits + 2), beta_Al);
//            mlmap_evaluation_N(bAs, 4, num_terms >> (logN + 1 + log_bits + 2), j >> (logN + 1 + log_bits + 2), beta_As);
//            mlmap_evaluation_N(bBl, 4, 0, (j >> (log_bits + 2)) & (N * 4 - 1), beta_Bl);
//            mlmap_evaluation_N(bBs, 4, 0, (j >> (log_bits + 2)) & (N * 2 - 1), beta_Bs);
//            mlmap_evaluation_N(bC, 4, 0, j & (8 * ubits - 1), beta_C);
//
//            mlmap_evaluation_N(c10, 4, num_terms, j, V_1c0);
//            mlmap_evaluation_N(c11, 4, num_terms, j, V_1c1);
//            mlmap_evaluation_N(c00, 4, num_terms, j, V_0c0);
//            mlmap_evaluation_N(c01, 4, num_terms, j, V_0c1);
//
//            mlmap_evaluation_N(a1, 3, 0, j & (8 * ubits - 1), alpha1);
//            mlmap_evaluation_N(a2, 3, 0, j & (8 * ubits - 1), alpha2);
//            mlmap_evaluation_N(a3, 3, 0, j & (8 * ubits - 1), alpha3);
//
//            mlmap_evaluation_N(tl0, 3, 0, (j >> (log_bits + 2)) & (N * 4 - 1), tau_l0);
//            mlmap_evaluation_N(tl1, 3, 0, (j >> (log_bits + 2)) & (N * 4 - 1), tau_l1);
//            mlmap_evaluation_N(ts, 3, 0, (j >> (log_bits + 2)) & (N * 2 - 1), tau_s);
//
//
//            //update poly
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, tl0[k]);
//                mod_mult(tmp1, tmp1, c00[k]);
//                mod_add(poly_0d[k], poly_0d[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, tl1[k]);
//                mod_mult(tmp1, tmp1, c00[k]);
//                mod_add(poly_0s0[k], poly_0s0[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, ts[k]);
//                mod_mult(tmp1, tmp1, c01[k]);
//                mod_add(poly_0s1[k], poly_0s1[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_sub(tmp1, a2[k], a1[k]);
//                mod_add(tmp1, tmp1, P);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, a3[k]);
//                mod_mult(tmp1, tmp1, ts[k]);
//                mod_mult(tmp1, tmp1, c01[k]);
//                mod_add(poly_0r[k], poly_0r[k], tmp1);
//            }
//            for (int k = 0; k < 4; k ++) {
//                mod_sub_ui(tmp1, c00[k], 1);
//                mod_mult(tmp1, tmp1, c00[k]);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, tl[k]);
//                mod_add(poly_0b0[k], poly_0b0[k], tmp1);
//            }
//            for (int k = 0; k < 4; k ++) {
//                mod_sub_ui(tmp1, c01[k], 1);
//                mod_mult(tmp1, tmp1, c01[k]);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, ts[k]);
//                mod_add(poly_0b1[k], poly_0b1[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, tl0[k]);
//                mod_mult(tmp1, tmp1, c10[k]);
//                mod_add(poly_1d[k], poly_1d[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, tl1[k]);
//                mod_mult(tmp1, tmp1, c10[k]);
//                mod_add(poly_1s0[k], poly_1s0[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_mult(tmp1, bAl[k], a3[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, ts[k]);
//                mod_mult(tmp1, tmp1, c11[k]);
//                mod_add(poly_1s1[k], poly_1s1[k], tmp1);
//            }
//            for (int k = 0; k < 3; k ++) {
//                mod_sub(tmp1, a2[k], a1[k]);
//                mod_add(tmp1, tmp1, P);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_mult(tmp1, tmp1, a3[k]);
//                mod_mult(tmp1, tmp1, ts[k]);
//                mod_mult(tmp1, tmp1, c11[k]);
//                mod_add(poly_1r[k], poly_1r[k], tmp1);
//            }
//            for (int k = 0; k < 4; k ++) {
//                mod_sub_ui(tmp1, c10[k], 1);
//                mod_mult(tmp1, tmp1, c10[k]);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBl[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_add(poly_1b0[k], poly_1b0[k], tmp1);
//            }
//            for (int k = 0; k < 4; k ++) {
//                mod_sub_ui(tmp1, c11[k], 1);
//                mod_mult(tmp1, tmp1, c11[k]);
//                mod_mult(tmp1, tmp1, bAl[k]);
//                mod_mult(tmp1, tmp1, bBs[k]);
//                mod_mult(tmp1, tmp1, bC[k]);
//                mod_add(poly_1b1[k], poly_1b1[k], tmp1);
//            }
// 
//        }
//
//        num_terms >>= 1;
//        update_V(beta_Al, )
//        update_V(beta_As, )
//        update_V(V_0c0, )
//        update_V(V_0c1, )
//        update_V(V_1c0, )
//        update_V(V_1c1, )
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
//    stat = stat ? sum_check_verification(tmp2, poly_r, vrr, 3, 1 + log_num + logN + log_bits, r_c, "Vr") : 0;
//    mod_sub(tmp1, a1r, a2r);
//    mod_mult(tmp1, tmp1, bhr);
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
//    stat = stat ? sum_check_verification(tmp2, poly_d, vdr, 3, 1 + log_num + logN + log_bits, r_c, "Vd") : 0;
//    mod_mult(tmp1, tr, bhr);
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
//    return stat;
//}
//


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
