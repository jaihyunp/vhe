/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *Currently working: this script may not work
 *vc7.c is the latest stable script
 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*************************************
//Jai Hyun Park
//October 7, 2019.
//Implementation of fast rounding circuit using commitments
**************************************/
#include <gmp.h>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>


static mpz_t PRIME;
typedef unsigned long long uint64;


void mod(mpz_t rop, const mpz_t x)
{
    mpz_mod(rop, x, PRIME);
}

void mod_add(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_add(rop, x, y);
    mod(rop, rop);
}
void mod_sub(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_sub(rop, x, y);
    mod(rop, rop);
}
void mod_mult(mpz_t rop, const mpz_t x, const mpz_t y)
{
    mpz_mul(rop, x, y);
    mod(rop, rop);
}
inline void mod_1neg(mpz_t rop, const mpz_t r)
{
    mpz_ui_sub(rop, 1, r);
    mod(rop, rop);
}
void mod_pow(mpz_t rop, const mpz_t x, uint64 b)
{
    uint64 d = 1;
    mpz_t base, res;
    mpz_inits(base, res, NULL);

    mpz_set_ui(res, 1);
    mpz_set(base, x);

    while(d <= b) {
        if(d & b)
            mod_mult(res, res, base);
        d = d << 1;
        if(d <= b)
            mod_mult(base, base, base);
    }

    mod(rop, res);
    mpz_clears(base, res, NULL);
}

void euclidean_algo(const mpz_t u, mpz_t u1, mpz_t u2, mpz_t u3)
{
    mpz_t v1, v2, v3, q, t1, t2, t3, tmp;
    mpz_inits(v1, v2, v3, q, t1, t2, t3, tmp, NULL);

    mpz_set_ui(u1, 1);
    mpz_set_ui(u2, 0);
    mpz_set(u3, u);

    mpz_set_ui(v1, 0);
    mpz_set_ui(v2, 1);
    mpz_set(v3, PRIME);

    do {
        mpz_fdiv_q(q, u3, v3);

        mod_mult(tmp, q, v1);
        mod_sub(t1, u1, tmp);

        mod_mult(tmp, q, v2);
        mod_sub(t2, u2, tmp);

        mod_mult(tmp, q, v3);
        mod_sub(t3, u3, tmp);
        
        mpz_set(u1, v1);
        mpz_set(u2, v2);
        mpz_set(u3, v3);
        
        mpz_set(v1, t1);
        mpz_set(v2, t2);
        mpz_set(v3, t3);
    } while(mpz_sgn(v3) != 0);

    mpz_clears(v1, v2, v3, q, t1, t2, t3, tmp, NULL);
}

void mod_inv(mpz_t rop, const mpz_t a)
{
    mpz_t u1, u2, u3;
    mpz_inits(u1, u2, u3, NULL);

    euclidean_algo(a, u1, u2, u3);

    if(!mpz_cmp_ui(u3, 1))
        mod(rop, u1);
    else
        mpz_set_ui(rop, 0);

    mpz_clears(u1, u2, u3, NULL);
}


void eval_poly(mpz_t rop, const mpz_t *coeff, const int degree, const mpz_t val)
{
    mpz_t base, tmp;
    mpz_inits(base, tmp, NULL);

    mpz_set(base, val);
    mpz_set(rop, coeff[0]);

    for(int i = 1; i <= degree; i ++) {
        mod_mult(tmp, base, coeff[i]);
        mod_add(rop, rop, tmp);
        if(i < degree)
            mod_mult(base, base, val);
    }

    mpz_clears(base, tmp , NULL);
    return;
}

void eval_polys(mpz_t *evs, const mpz_t **coeff, const int n, const int degree, const mpz_t val)
{
    for(int i = 0; i < n; i ++)
        eval_poly(evs[i], coeff[i], degree, val);
    return;
}



void extrapolate(mpz_t rop, const mpz_t *vec, const int n, const mpz_t r)
{
    mpz_t mult, tmp, inv;
    mpz_inits(mult, tmp, inv, NULL);

    mpz_set_ui(rop, 0);
    for(int i = 0; i < n; i ++) {
        mpz_set_ui(mult, 1);
        for(int j = 0; j < n; j ++) {
            if(i != j){
                mpz_set_si(tmp, i - j);
                mod_inv(inv, tmp);
                mpz_sub_ui(tmp, r, (unsigned int) j);
                mod_mult(mult, mult, tmp);
                mod_mult(mult, mult, inv);
            }
        }
        mod_mult(tmp, mult, vec[i]);
        mod_add(rop, rop, tmp);
    }

    mpz_clears(mult, tmp, inv, NULL);
}


void update_V(mpz_t *V, const int num_new, const mpz_t ri)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    for(int i = 0; i < num_new; i ++){
        mod_1neg(tmp1, ri);
        mod_mult(tmp1, tmp1, V[i]);

        mod_mult(tmp2, ri, V[i + num_new]);

        mod_add(V[i], tmp1, tmp2);
    }
    mpz_clears(tmp1, tmp2, NULL);
}

void evaluate_V(mpz_t rop, const mpz_t *V_in, const int d, const mpz_t *r)
{
    uint64 n = 1ULL << d;
    mpz_t *V = (mpz_t*) malloc(sizeof(mpz_t) << d);
    for(uint64 i = 0; i < n; i ++)
        mpz_init_set(V[i], V_in[i]);

    for(int round = 0; round < d; round ++) {
        update_V(V, n >> (round + 1), r[d - round - 1]);
    }

    mod(rop, V[0]);
    for(uint64 i = 0; i < n; i ++)
        mpz_clear(V[i]);
    free(V);
}

void initialize_beta(mpz_t* betavals, const int d, const mpz_t* z)
{
    mpz_t zval, oldval, tmp;
    uint64 two_to_k = 1;

    mpz_inits(zval, oldval, tmp, NULL);

    mpz_set_ui(betavals[0], 1);
    for(int k = 0; k < d; k ++) {
        mpz_set(zval, z[k]);
        for(uint64 i = 0; i < two_to_k; i ++) {
            mpz_set(oldval, betavals[i]);
            mod_1neg(tmp, zval);
            mod_mult(betavals[i], oldval, tmp);
            mod_mult(betavals[i + two_to_k], oldval, zval);
        }
        two_to_k <<= 1;
    }

    mpz_clears(zval, oldval, tmp, NULL);

}
void evaluate_beta(mpz_t rop, const mpz_t* z, const mpz_t* r, const int d)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    mpz_set_ui(rop, 1);
    for(int i = 0; i < d; i ++){
        mod_1neg(tmp1, z[i]);
        mod_1neg(tmp2, r[i]);
        mod_mult(tmp1, tmp1, tmp2);//tmp1<-(1-zi)*(1-ri)

        mod_mult(tmp2, z[i], r[i]);//tmp2<-zi*ri

        mod_add(tmp1, tmp1, tmp2);//tmp1<-(1-zi)*(1-ri)+zi*ri

        mod_mult(rop, rop, tmp1);
    }
    mpz_clears(tmp1, tmp2, NULL);
}

#define ERROR 0x0
int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);


    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check %s %s\n", str, mpz_get_str(NULL, 10, tmp), mpz_get_str(NULL, 10, ri));
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

    if(stat)
        mpz_set(rop, extrap_val);
    else
        mpz_set_ui(rop, ERROR);

    mpz_clears(tmp, extrap_val, NULL);

    return stat;
}

void linear_evaluation(mpz_t *v, const int n, const uint64 num_terms, const int pivot, const mpz_t* V)
{
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    mpz_set(v[0], V[pivot]);
    mpz_set(v[1], V[pivot + (num_terms >> 1)]);

    for(int i = 2; i < n; i ++) {
        mpz_mul_si(tmp1, v[1], i);
        mpz_mul_si(tmp2, v[0], 1 - i);
        mod_add(v[i], tmp1, tmp2);
    }

    mpz_clears(tmp1, tmp2, NULL);
}


void mod_inv_pow2(mpz_t rop, const int k)
{
    mpz_t tmp;
    mpz_init_set_ui(tmp, 2);
    mod_pow(tmp, tmp, k);
    mod_inv(rop, tmp);
    mpz_clear(tmp);
}

void evaluate_alpha(mpz_t rop, const uint64 e_in, const mpz_t* r, const int d)
{
    mpz_t l, u, rval, tmp, tmp1, tmp2;
    uint64 e = e_in - 1;

    if((e >= (1ULL << d)) || (e_in == 0)) {
        mpz_set_ui(rop, 0);
        return;
    }
    
    mpz_inits(l, u, rval, tmp, tmp1, tmp2, NULL);
    mpz_set_ui(l, 1);
    mpz_set_ui(u, 0);
    for(int i = d - 1; i >= 0; i --) {
        mpz_set_ui(tmp1, 2);
        mod_pow(tmp1, tmp1, 1ULL << i);
        mod_mult(rval, r[i], tmp1); //rval<-r[i]*2^(2^i)


        if((e >> i) & 1) {
            mod_1neg(tmp1, r[i]);
            mod_mult(tmp1, tmp1, l);
            mod_mult(tmp2, rval, l);
            mod_add(l, tmp1, tmp2);

            mod_mult(tmp2, rval, u);
            mod_add(u, tmp1, tmp2);
        } else {
            mod_1neg(tmp1, r[i]);
            mod_mult(tmp1, tmp1, l);
            mod_mult(tmp2, rval, u);
            mod_add(l, tmp1, tmp2);

            mod_1neg(tmp1, r[i]);
            mod_mult(tmp1, tmp1, u);
            mod_add(u, tmp1, tmp2);
        }
    }

    mpz_set(rop, l);
    mpz_clears(l, u, rval, tmp1, tmp2, NULL);
}



void evaluate_tau(mpz_t rop, const mpz_t val, const mpz_t *r, const int d)
{
    mpz_t base, tmp, one;
    mpz_set_ui(rop, 1);
    mpz_init_set(base, val);
    mpz_init(tmp);
    mpz_init_set_ui(one, 1);
    for(int i = 0; i < d; i ++){
        mod_sub(tmp, base, one);
        mod_mult(tmp, tmp, r[d - 1 - i]);
        mod_add(tmp, tmp, one);
        mod_mult(rop, rop, tmp);
        if(i < (d - 1))
            mod_mult(base, base, base);
    }
    mpz_clears(base, tmp, one, NULL);
}
void initialize_alpha(mpz_t *alpha, const int l, const int e)
{
    mpz_t base, two;
    mpz_init_set_ui(base, 1);
    mpz_init_set_ui(two, 2);
    for(int i = 0; i < (1 << l); i ++){
        if(i < e) {
            mpz_set(alpha[i], base);
            mod_mult(base, base, two);
        } else {
            mpz_set_ui(alpha[i], 0);
        }
    }
    mpz_clears(base, two, NULL);
}
void initialize_tau(mpz_t *tau, const int l, const mpz_t val)
{
    mpz_t base;
    mpz_init_set_ui(base, 1);
    for(int i = 0; i < (1 << l); i ++){
        mpz_set(tau[i], base);
        mod_mult(base, base, val);
    }
    mpz_clear(base);
}
// log # poly = ln
// log degree = ld
// d1 bits -> d2 bits
int sum_check_poly_rounding(mpz_t rop, mpz_t *C, const mpz_t ri, const int ln, const int ld, const int d1, const int d2, const mpz_t *r_n, const mpz_t *r_d, const mpz_t *r_b, const mpz_t val, const mpz_t *z_n, const mpz_t *z_db)
{
    int stat = 1;
    int lb = 0;
    while((d1 - 1) >> ++ lb);
    int l = ln + ld + lb;
    mpz_t tmp1, tmp2, cr, a1r, a2r, b1r, b2r, tr;
    mpz_inits(tmp1, tmp2, cr, a1r, a2r, b1r, b2r, tr, NULL);
    mpz_t *r = (mpz_t*) malloc(sizeof(mpz_t) * (ln + ld + lb)),
          *z = (mpz_t*) malloc(sizeof(mpz_t) * (ln + ld + lb));
    for(int i = 0; i < ln; i ++)
        mpz_init_set(r[i], r_n[i]);
    for(int i = 0; i < ld; i ++)
        mpz_init_set(r[i + ln], r_d[i]);
    for(int i = 0; i < lb; i ++)
        mpz_init_set(r[i + ln + ld], r_b[i]);
    for(int i = 0; i < ln; i ++)
        mpz_init_set(z[i], z_n[i]);
    for(int i = 0; i < ld + lb; i ++)
        mpz_init_set(z[i + ln], z_db[i]);
    
    //commit
    evaluate_V(cr, C, l, r);
    evaluate_alpha(a1r, d1, r_b, lb);
    evaluate_alpha(a2r, d1 - d2, r_b, lb);
    evaluate_beta(b1r, z, r_n, ln);
    evaluate_beta(b2r, z, r, l);
    evaluate_tau(tr, val, r_d, ld); 

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
    initialize_alpha(alpha1, lb, d1);
    initialize_alpha(alpha2, lb, d1 - d2);
    initialize_tau(tau, ld, val);

    mpz_t **poly_b = (mpz_t**) malloc(sizeof(mpz_t*) * l),
          **poly_d = (mpz_t**) malloc(sizeof(mpz_t*) * l),
          **poly_r = (mpz_t**) malloc(sizeof(mpz_t*) * l);
    for(int i = 0; i < l; i ++){
        poly_b[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        poly_d[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        poly_r[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        for(int j = 0; j < 4; i ++){
            mpz_init_set_ui(poly_b[i][j], 0);
        }
        for(int j = 0; j < 3; i ++){
            mpz_init_set_ui(poly_d[i][j], 0);
            mpz_init_set_ui(poly_r[i][j], 0);
        }
    }

    //proof
    uint64 num_terms = 1ULL << l;
    for(int round = 0; round < l; round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++) {

            linear_evaluation(b1, 3, num_terms >> (ld + lb), j >> (ld + lb), beta1);
            linear_evaluation(t, 3, (num_terms >> lb) & ((1ULL << ld) - 1), (j >> lb) & ((1ULL << ld) - 1), tau);
            linear_evaluation(a1, 3, num_terms & ((1ULL << lb) - 1), j & ((1ULL << lb) - 1), alpha1);
            linear_evaluation(a2, 3, num_terms & ((1ULL << lb) - 1), j & ((1ULL << lb) - 1), alpha2);
            
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
            update_V(tau, num_terms >> ld, r[l - 1 - round]);
        } else {
            update_V(alpha1, num_terms, r[l - 1 - round]);
            update_V(alpha2, num_terms, r[l - 1 - round]);
        }
        update_V(C, num_terms, r[l - 1 - round]);
        update_V(beta2, num_terms, r[l - 1 - round]);
        
    }


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


    //return
    if(stat)
        mod_add(rop, poly_d[0][0], poly_d[0][1]);
    else
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
        mpz_init(beta1[i]);
    for(uint64 i = 0; i < (1ULL << l); i ++)
        mpz_init(beta2[i]);
    for(uint64 i = 0; i < (1ULL << ld); i ++)
        mpz_init(tau[i]);
    for(uint64 i = 0; i < (1ULL << lb); i ++)
        mpz_inits(alpha1[i], alpha2[i], NULL);
    free(beta1);
    free(beta2);
    free(alpha1);
    free(alpha2);
    free(tau);

    for(int i = 0; i < l; i ++){
        poly_b[i] = (mpz_t*) malloc(sizeof(mpz_t) * 4);
        poly_d[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        poly_r[i] = (mpz_t*) malloc(sizeof(mpz_t) * 3);
        for(int j = 0; j < 4; i ++){
            mpz_init_set_ui(poly_b[i][j], 0);
        }
        for(int j = 0; j < 3; i ++){
            mpz_init_set_ui(poly_d[i][j], 0);
            mpz_init_set_ui(poly_r[i][j], 0);
        }
    }

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

int main(int argc, char** argv)
{

    if((argc != 4)) {
        std::cout << "argc!= 2. Command line arg should be log(n) for nxn matrix multiplication\n";
        exit(1);
    } 
  
    clock_t t = clock();
    int stat = 0;
    int d = atoi(argv[1]);
    int n = 1 << d;
    int num_digits = atoi(argv[2]), num_rdigits = atoi(argv[3]), del = num_digits - num_rdigits;
    int l = 0;
    while((num_digits - 1) >> ++ l);
    mpz_t tmp1, ri, ri0;
    mpz_inits(tmp1, ri, ri0, NULL);

//    mpz_init_set_str(PRIME, "1FFFFFFFFFFFFFFF", 16);
    mpz_init(PRIME);
//    mpz_ui_pow_ui(PRIME, 2, 521);
//    mpz_ui_pow_ui(PRIME, 2, 61);
//    mpz_sub_ui(PRIME, PRIME, 1);
//    mpz_set_str(PRIME, "114993957235856291048434975234089602617",10 );
    mpz_set_str(PRIME, "205311158772819412817284986118022354041", 10);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));
  
    mpz_t* z = (mpz_t*) malloc((d + l)* sizeof(mpz_t));
    mpz_t* r = (mpz_t*) malloc((d + l) * sizeof(mpz_t));
    mpz_t* r_tail = (mpz_t*) malloc(d * sizeof(mpz_t));

    for(int i = 0; i < d + l; i ++){
        mpz_init_set_ui(r[i], rand());
        mod(r[i], r[i]);
//        printf("%d: %s\n", i, mpz_get_str(NULL, 16, r[i]));
//        mpz_init_set_ui(r[i], 100 + i);
    }
    for(int i = 0; i < d + l; i ++){
//        mpz_init_set_ui(z[i], 120 + i);
        mpz_init_set_ui(z[i], rand());
        mod(z[i], z[i]);
//        printf("%d: %s\n", i, mpz_get_str(NULL, 16, z[i]));
    }
    for(int i = 0; i < d; i ++) {
        mpz_init_set(r_tail[i], r[i + l]);
    }

    //run through entire Muggles protocol with prover
    mpz_t* V1 = (mpz_t*) malloc(n * sizeof(mpz_t));
    mpz_t* V0 = (mpz_t*) malloc(n * sizeof(mpz_t));
    mpz_t* C = (mpz_t*) malloc(sizeof(mpz_t) << (d + l));    
    for(int i = 0; i < n; i ++) {
        mpz_init_set_ui(V1[0], 0);
        mpz_init_set_ui(V0[0], 0);
    }
    for(int i = 0; i < (1 << (d + l)); i ++) {
        mpz_init_set_ui(C[i], 0);
    }
    
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < num_digits; j ++) {
            unsigned int random_bit = rand() & 1;
            
            if(random_bit) {
                mpz_set_ui(C[(i << l) + j], 1);

                mpz_set_ui(tmp1, 2);
                mod_pow(tmp1, tmp1, j);
                mpz_add(V1[i], V1[i], tmp1);

                if(j >= del) {
                    mpz_set_ui(tmp1, 2);
                    mod_pow(tmp1, tmp1, j - del);
                    mpz_add(V0[i], V0[i], tmp1);
                }
            }
            
        }
    }

    printf("%f sec til main init\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
  
    evaluate_V(ri, V0, d, r_tail);
    evaluate_V(ri0, V1, d, r_tail);
    stat = sum_check_rounding_fast(ri, C, ri, d, num_digits, num_rdigits, r, z);
    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
    } else {
        printf("Failed inside the layer\n");
    }
    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
  
    return 1;
}


