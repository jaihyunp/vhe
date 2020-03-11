#include "mlmap.h"

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
        mod_pow_ui(tmp1, tmp1, 1ULL << i);
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

void evaluate_tau2(mpz_t rop, const mpz_t val, const mpz_t *r, const int d)
{
    mpz_t base, tmp, one;
    mpz_set_ui(rop, 1);
    mpz_init_set(base, val);
    mpz_init(tmp);
    mpz_init_set_ui(one, 1);

    mod_sub(rop, one, r[d - 1]);
    mod_sub(rop, rop, r[d - 1]);
    for(int i = 0; i < d - 1; i ++) {
        mod_sub(tmp, base, one);
        mod_mult(tmp, tmp, r[i]);
        mod_add(tmp, tmp, one);
        mod_mult(rop, rop, tmp);
        if(i < (d - 2))
            mod_mult(base, base, base);
 
    }
    mpz_clears(base, tmp, one, NULL);
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
        mod_mult(tmp, tmp, r[i]);
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

void initialize_tau2(mpz_t *tau, const int l, const mpz_t val)
{
    mpz_set_ui(tau[0], 1);
    mpz_sub_ui(tau[1 << (l - 1)], PRIME, 1);
//    mpz_set_si(tau[1 << (l - 1)], -1);
//    mod(tau[1 << (l - 1)], tau[1 << (l - 1)]);

    for(int i = 1; i < (1 << (l - 1)); i ++){
        mod_mult(tau[i], tau[i - 1], val);
        mod_mult(tau[i + (1 << (l - 1))], tau[i - 1 + (1 << (l - 1))], val);
    }
}
void initialize_tau(mpz_t *tau, const int l, const mpz_t val)
{
    mpz_set_ui(tau[0], 1);
    for(int i = 1; i < (1 << l); i ++)
        mod_mult(tau[i], tau[i - 1], val);
}

