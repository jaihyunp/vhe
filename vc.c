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

//#define MASK 4294967295 //2^32-1
//#define PRIME 2305843009213693951 //2^61-1
//#define PROOF 1

static mpz_t PRIME;

typedef unsigned long long uint64;

void mod(mpz_t rop, const mpz_t x)
{
    mpz_mod(rop, x, PRIME);
}
/*
//and if so subtract p, but I don't want to pay that
//efficiency hit here, so the user should just be aware of this.
inline uint64 myMod(uint64 x)
{
    uint64 res = (x >> 61) + (x & PRIME);
    if(res < PRIME)
        return res;
    else
        return res - PRIME;
}
*/


/*
//computes x^b
inline uint64 myPow(uint64 x, uint64 b)
{
    uint64 res = 1, base = x, d = 1;//base = x ** d, d = 2 ** i
    while(d <= b) {
        if(d & b)
            res *= base;
        base *= base;
        d = d << 1;
    }
    return res;
}
*/

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
/*
//efficient modular multiplication function mod 2^61-1
inline uint64 myModMult(uint64 x, uint64 y)
{
	uint64 hi_x = x >> 32; 
	uint64 hi_y = y >> 32;
	uint64 low_x = x & MASK;
	uint64 low_y = y & MASK;

	//since myMod might return something slightly large than 2^61-1,  
	//we need to multiply by 8 in two pieces to avoid overflow.
	uint64 piece1 = myMod((hi_x * hi_y) << 3);
	uint64 z = (hi_x * low_y + hi_y * low_x);
	uint64 hi_z = z >> 32;
	uint64 low_z = z & MASK;

	//Note 2^64 mod (2^61-1) is 8
	uint64 piece2 = myMod((hi_z << 3) + myMod((low_z << 32)));
	uint64 piece3 = myMod(low_x * low_y);
	uint64 result = myMod(piece1 + piece2 + piece3);
	
    return result;
}
*/

/*
//computes b^e mod p using repeated squampz_t. p should be 2^61-1
inline uint64 myModPow(uint64 b, uint64 e)
{
    uint64 result;
    if(e == 1)
        return b;
    if(e == 0)
        return 1;
    if((e & 1) == 0) {
        result = myModPow(b, (e >> 1));
        return myModMult(result, result);
    } else {
        return myModMult(myModPow(b, e - 1), b);
    }
}
*/


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
/*
//Performs Extended Euclidean Algorithm
//Used for computing multiplicative inverses mod p
void extEuclideanAlg(uint64 u, uint64* u1, uint64* u2, uint64* u3)
{
    *u1 = 1;
    *u2 = 0;
    *u3 = u;
    uint64 v1 = 0;
    uint64 v2 = 1;
    uint64 v3 = PRIME;
    uint64 q;
    uint64 t1;
    uint64 t2;
    uint64 t3;
    do
    {
        q = *u3 / v3;
        //t1 = *u1 + p - q * v1;
        //t2 = *u2 + p - q*v2;
        //t3 = *u3 + p - q*v3;
        t1 = myMod((*u1) + PRIME - myModMult(q, v1));
        t2 = myMod((*u2) + PRIME - myModMult(q, v2));
        t3 = myMod((*u3) + PRIME - myModMult(q, v3));
        (*u1) = v1;
        (*u2) = v2;
        (*u3) = v3;
        v1 = t1;
        v2 = t2;
        v3 = t3;
    } while(v3 != 0);
}
*/


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
/*
//Computes the modular multiplicative inverse of a modulo m,
//using the extended Euclidean algorithm
//only works for p=2^61-1
inline uint64 inv(uint64 a)
{
    uint64 u1, u2, u3, res;
    extEuclideanAlg(a, &u1, &u2, &u3);
    if(u3 ==  1)
          return myMod(u1);
    else
          return 0;
}
*/

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
/*
//extrapolate the polynomial implied by vector vec of length n to location r
uint64 extrap(uint64* vec, uint64 n, uint64 r)
{
    uint64 result = 0;
    uint64 mult = 1;
    for(uint64 i = 0; i < n; i ++) {
        mult = 1;
        for(uint64 j = 0; j < n; j ++) {
            if (i > j) {
                if(r > j)
                    mult = myModMult(myModMult(mult, r - j), inv(i - j));
                else
                    mult = myModMult(myModMult(mult, PRIME + r - j), inv(i - j));
            }
            if (i < j) {
                if(r > j)
                    mult = myModMult(myModMult(mult, r - j), inv(PRIME + i - j));
                else
                    mult = myModMult(myModMult(mult, PRIME + r - j), inv(PRIME + i - j));
            }
        }
        result = myMod(result + myModMult(mult, vec[i]));
    }
    return result;
}
*/


inline void mod_1neg(mpz_t rop, const mpz_t r)
{
    mpz_ui_sub(rop, 1, r);
    mod(rop, rop);
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
/*
inline void updateV(uint64* V, int num_new, uint64 ri)
{
    for(int i = 0; i < num_new; i ++) 
		V[i] = myMod(myModMult(V[i], 1 + PRIME - ri) + myModMult(V[i + num_new], ri));
}
*/


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
/*
uint64 evaluate_V(uint64* V_in, int d, uint64 *r)
{
    uint64 n = 1ULL << d, ans;
    uint64 *V = (uint64*) malloc(sizeof(uint64) << d);

    for(uint64 i = 0; i < n; i ++)
        V[i] = V_in[i];

    for(int round = 0; round < d; round ++) {
        updateV(V, n >> (round + 1), r[d - round - 1]);
    }
    ans = V[0];

    free(V);
    return ans;
}
*/


void initialize_betavals(mpz_t* betavals, const int d, const mpz_t* z)
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
/*
//sets betavals(p) = \prod_{i=0}^{d-1} (zipi + (1-zi)(1-pi))
void initialize_betavals(uint64* betavals, int d, uint64* z)
{
	uint64 zval, oldval;
	int two_to_k = 1;
	
    betavals[0] = 1;
	for(int k = 0; k < d; k ++) {
		zval = z[k];
		for(int i = 0; i < two_to_k; i ++) {
			oldval = betavals[i];
			betavals[i] = myModMult(oldval, 1 + PRIME - zval);
			betavals[i + two_to_k] = myModMult(oldval, zval);
		}
		two_to_k <<= 1;
	}
}
*/


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
/*
inline uint64 evaluate_beta(uint64 *z, uint64 *r, int d)
{
    uint64 res = 1;
    for(int i = 0; i < d; i ++)
        res = myModMult(res, myMod(myModMult(1 + PRIME - z[i], 1 + PRIME - r[i]) + myModMult(z[i], r[i])));
    return res;
}
*/

#define ERROR 0x0
int sum_check_verification(mpz_t rop, mpz_t** poly, const mpz_t ri, const int degree, const int num_rounds, const mpz_t *r, const char *str = NULL)
{
    int stat = 1;
    mpz_t tmp, extrap_val;
    mpz_inits(tmp, extrap_val, NULL);


    mod_add(tmp, poly[0][0], poly[0][1]);
    if(mpz_cmp(tmp, ri)) {
        printf("Fail Sumcheck %s :: 1st sum check\n", str);
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
/*
inline uint64 sum_check_verification(uint64 **poly, uint64 ri, int degree, int num_rounds, uint64 *r, const char *str = NULL)
{
    if(myMod(poly[0][0] + poly[0][1]) != ri){
        printf("Fail sumcheck %s :: 1st sum check %lld %lld\n", str, myMod(poly[0][0] + poly[0][1]), ri);
        return -1;
    }

    uint64 extrap_val = extrap(poly[0], degree, r[num_rounds - 1]);
    for(int round = 1; round < num_rounds; round ++){
        if(myMod(poly[round][0] + poly[round][1]) != extrap_val){
            printf("Fail sumcheck %s :: %dth round\n", str, round);
            return -1;
        } else {
            extrap_val = extrap(poly[round], degree, r[num_rounds - 1 - round]);
        }
    }

    return extrap_val;
}
*/


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
/*
//evaluate V(0, p), V(1, p), .., V(n, p)
inline void linear_evaluation(uint64* v, int n, uint64 num_terms, int pivot, uint64* V)
{
    v[0] = V[pivot];
    v[1] = V[pivot + (num_terms >> 1)];
    if(n < 8){
        for(int i = 2; i < n; i ++){
            v[i] = myMod(i * v[1] + (i - 1) * (PRIME - v[0]));
        }
    } else {
        for(int i = 2; i < n; i ++)
            v[i] = myMod(myModMult(i, v[1]) + myModMult(1 + PRIME - i, v[0]));
    }
}
*/

/*
//old version
uint64 sum_check_rounding(uint64* V, uint64 **C, uint64 ri, int d, int num_digits, int num_rdigits, uint64 *r, uint64 *z)
{
    uint64 ***poly, *betavals;
    uint64 stat = 0, res = 0;
    clock_t t = clock();

    //Assumption: Commitment enables us to evaluate C[i] with verification
    uint64 *cr = (uint64*) malloc(num_digits * sizeof(uint64));
    uint64 *cz = (uint64*) malloc(num_digits * sizeof(uint64));
    for(int i = 0; i < num_digits; i ++){
        cr[i] = evaluate_V(C[i], d, r);
        cz[i] = evaluate_V(C[i], d, z);
    }

    printf("%f sec til cr\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();

    //initialize poly
    poly = (uint64***) malloc(num_digits * sizeof(uint64**));
    for(int i = 0; i < num_digits; i ++){
        poly[i] = (uint64**) malloc(d * sizeof(uint64*));
        for(int j = 0; j < d; j ++){
            poly[i][j] = (uint64*) calloc(4, sizeof(uint64));
        }
    }

    //initialize betavals
    betavals = (uint64*) malloc((1ULL << d) * sizeof(uint64));
    //uint64 br = evaluate_V(betavals, d, r);
    uint64 br = evaluate_beta(z, r, d);

    printf("%f sec til br\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();


    initialize_betavals(betavals, d, z);
    //prover do:
    uint64 b[4], c[4];
    uint64 num_terms =  1ULL << d;
    for(int round = 0; round < d; round ++){

        for(int i = 0; i < num_digits; i ++){
            for(uint64 j = 0; j < num_terms / 2; j ++){

                linear_evaluation(b, 4, num_terms, j, betavals);
                linear_evaluation(c, 4, num_terms, j, C[i]);
                
                for(int k = 0; k < 4; k ++)
                    poly[i][round][k] = myMod(poly[i][round][k] + myModMult(myModMult(b[k], c[k]), c[k] + PRIME - 1));

            }

    		updateV(C[i], num_terms >> 1, r[d - round - 1]);
        }
            
        num_terms >>= 1;
		updateV(betavals, num_terms, r[d - round - 1]);
        
    }
    for(int i = 0; i < num_digits; i ++)
        res = myMod(res + myModMult(1ULL << i, C[i][0]));
 
    printf("%f sec til pro\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();


    //verifier do:
    uint64 vr = 0;
    for(int i = num_digits - num_rdigits; i < num_digits; i ++)
        vr = myMod(vr + myModMult(1ULL << (i - num_digits + num_rdigits), cz[i]));
    if(vr != ri){
        printf("Fail :: not match with the previous layer\n");
        stat = 1;
    }
    
    if(!stat) {
        for(int i = 0; i < num_digits; i ++){
            if((!stat) && (sum_check_verification(poly[i], 0, 4, d, r, "V_b") != myModMult(br, myModMult(cr[i], cr[i] + PRIME - 1)))){
                printf("Fail :: %dth digit is_bit sumcheck does not match with the commitment\n", i);
                stat = 1;
                break;
            }
        }
    }

    printf("%f sec til ver\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();

    //free
    for(int i = 0; i < num_digits; i ++){
        for(int j = 0; j < d; j ++){
            free(poly[i][j]);
        }
        free(poly[i]);
    }
    free(poly);
    free(betavals);
    free(cr);
    free(cz);

    //return
    return stat ? 0 : res;
}
*/


void mod_inv_pow2(mpz_t rop, const int k)
{
    mpz_t tmp;
    mpz_init_set_ui(tmp, 2);
    mod_pow(tmp, tmp, k);
    mod_inv(rop, tmp);
    mpz_clear(tmp);
}
/*
inline uint64 inv_pow2(int k)
{
    return myMod(1ULL << (61 - k)); 
}
*/


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
/*
uint64 evaluate_alpha(int e_in, uint64 *r, int d)
{
    uint64 l = 1, u = 0, rval, e = e_in - 1, tmp;

    if((e >= (1ULL << d)) || (e_in == 0))
        return 0;

    for(int i = d - 1; i >= 0; i --){
        
        rval = myModMult(r[i], 1ULL << (1ULL << i));
        if((e >> i) & 1){
            tmp = myMod(myModMult(1 + PRIME - r[i], l) + myModMult(rval, l));
            u = myMod(myModMult(1 + PRIME - r[i], l) + myModMult(rval, u));
            l = tmp;
        } else {
            tmp = myMod(myModMult(1 + PRIME - r[i], l) + myModMult(rval, u));
            u = myMod(myModMult(1 + PRIME - r[i], u) + myModMult(rval, u));
            l = tmp;
        }
    }

    return l;
}
*/



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
    initialize_betavals(betavals, l + d, z); 
    
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

/*
//N = 2^d elements
//Rounding d1 digits to d2 digits
//Updates V_r, V, C, at the same time til z
uint64 sum_check_rounding_fast(uint64* V, uint64* C, uint64 ri, int d, int d1, int d2, uint64* r, uint64* z)
{
    int stat = 0;

    //d1, d2 are l-bits
    //For example, if a element from V is 60 bits and an element from V_r is 3 bits,
    //d1 = 60, d2 = 3 and l = 5; since 60 is 5 bits
    int l = 0;
    while((d1 - 1) >> ++ l);
    clock_t t = clock();

    //I omit commitment phase here
    //I assume there is a commitment algo fit our scheme
    //So C[r] is always available
    uint64 cr = evaluate_V(C, d + l, r);
    
    printf("%f sec til cr\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
    
    //construct and pre-compute the alpha (prover's task)
    uint64 *r_head = (uint64*) malloc(l * sizeof(uint64)),
           *r_tail = (uint64*) malloc(d * sizeof(uint64));
    for(int i = 0; i < l; i ++)
        r_head[i] = r[i];
    for(int i = 0; i < d; i ++)
        r_tail[i] = r[i + l];

    uint64 *alpha1 = (uint64*) calloc(1 << l, sizeof(uint64)),
           *alpha2 = (uint64*) calloc(1 << l, sizeof(uint64)),
           *betavals = (uint64*) malloc(sizeof(uint64) << (l + d));
    
    for(int i = 0; i < d1; i ++)
        alpha1[i] = myMod(1ULL << i);
    for(int i = d1 - d2; i < d1; i ++)
        alpha2[i] = myMod(1ULL << (i - d1 + d2));

    uint64 a1r = evaluate_alpha(d1, r_head, l);
    uint64 a2r = myModMult(a1r + PRIME - evaluate_alpha(d1 - d2, r_head, l), inv_pow2(d1 - d2));
    uint64 br = evaluate_beta(z, r, l + d);

    printf("%f sec til br\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
    
    initialize_betavals(betavals, l + d, z);
    
    
    //poly initialization
    uint64 **poly_b = (uint64**) malloc((d + l) * sizeof(uint64*)),
           **poly_d = (uint64**) malloc(l * sizeof(uint64*)),
           **poly_r = (uint64**) malloc(l * sizeof(uint64*));
    for(int i = 0; i < d + l; i ++)
        poly_b[i] = (uint64*) calloc(4, sizeof(uint64));
    for(int i = 0; i < l; i ++){
        poly_d[i] = (uint64*) calloc(4, sizeof(uint64));
        poly_r[i] = (uint64*) calloc(4, sizeof(uint64));
    }
    

    //prover do:
    //This scheme should be simultaneously executed with the previous step
    // So that it can use the same z, r (and hence the same betavals and ri)
    uint64 b[4], c[4];
    uint64 num_terms = 1 << (d + l);
    for(int round = 0; round < d; round ++) {
        for(uint64 j = 0; j < num_terms / 2; j ++){

            linear_evaluation(b, 4, num_terms, j, betavals);
            linear_evaluation(c, 4, num_terms, j, C);

            for(int k = 0; k < 4; k ++)
                poly_b[round][k] = myMod(poly_b[round][k] + myModMult(b[k], myModMult(c[k], c[k] + PRIME - 1)));

            //Co-execution with the previous step:
            // should update poly_vr's by using vr's
            // e.g.,
            //  if((num_terms & ((1 << l) - 1))) == 0){
            //      uint64 vr0 = Vr[j >> l], vr1 = [(j + tmp) >> l], .., vrn;
            //      for(int i = 0; i < n; i ++)
            //          poly_vr[round][i] += bi * poly(vri);
            //  }
        }

        num_terms >>= 1;
        updateV(betavals, num_terms, r_tail[d - 1 - round]);
        updateV(C, num_terms, r_tail[d - 1 - round]);

        //Co-execution with the previous step:
        // We should update Vr here
        // e.g.,
        //  updateV(Vr, num_terms >> l, r[d + l - 1 - round]);
    }

    printf("%f sec til pro1\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
 

    //verifier do:
    uint64 ri_b = sum_check_verification(poly_b, 0, 4, d, r_tail, "inter V_b");
    if(ri_b == 0){
        stat = 1;
        printf("Fail :: inter V_b\n");
    }

    printf("%f sec til ver1\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
    
    //Co-execution with the previous step:
    // Verifier investigates the process of the previous layer here.
    // When we co-execute, we verifies the previous layer and update the ri
    // e.g.,
    //  if(poly_vr[0][0] + poly_vr[0][1] - ri)
    //      return FALSE;
    //  extrap_val = extrap(poly_vr[0], n, r[d + l - 1]);
    //  for(int round = 1; round < d; round ++)
    //      if(poly_vr[round][0] + poly_vr[round][1] - extrap_val)
    //          return FALSE;
    //      else
    //          extrap_val = extrap(poly_vr[round], n, r[d + l - 1 - round]);
    //  ri = extrap_val;
    


    //prover do:
    for(int i = 0; i < l; i ++)
        for(int j = 0; j < 4; j ++)
            poly_b[i][j] = 0;

    if(!stat) {
        uint64 a1[3], a2[3];
        for(int round = 0; round < l; round ++) {
            for(uint64 j = 0; j < num_terms / 2; j ++) {
    
                linear_evaluation(b, 4, num_terms, j, betavals);
                linear_evaluation(c, 4, num_terms, j, C);
                linear_evaluation(a1, 3, num_terms, j, alpha1);
                linear_evaluation(a2, 3, num_terms, j, alpha2);

                for(int k = 0; k < 3; k ++){
                    poly_d[round][k] = myMod(poly_d[round][k] + myModMult(a1[k], c[k]));
                    poly_r[round][k] = myMod(poly_r[round][k] + myModMult(a2[k], c[k]));
                }
                for(int k = 0; k < 4; k ++){
                    poly_b[round][k] = myMod(poly_b[round][k] + myModMult(b[k], myModMult(c[k], c[k] + PRIME - 1)));
                }
            
            }

            num_terms >>= 1;
            updateV(betavals, num_terms, r_head[l - 1 - round]);
            updateV(C, num_terms, r_head[l - 1 - round]);
            updateV(alpha1, num_terms, r_head[l - 1 - round]);
            updateV(alpha2, num_terms, r_head[l - 1 - round]);
        }
    }

    printf("%f sec til pro2\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
 

    //verifier do:
    if ((!stat) && (myModMult(br, myModMult(cr, cr + PRIME - 1)) != sum_check_verification(poly_b, ri_b, 4, l, r_head, "V_b"))){
        stat = 1;
        printf("Fail :: V_b\n");
    }
    
    if ((!stat) && (myModMult(a2r, cr) != sum_check_verification(poly_r, ri, 3, l, r_head, "V_r"))){
        stat = 1;
        printf("Fail :: V_r\n");
    }

    stat = 0;
    uint64 res = myMod(poly_d[0][0] + poly_d[0][1]);
    if ((!stat) && (myModMult(a1r, cr) != sum_check_verification(poly_d, res, 3, l, r_head, "V_d"))){
        stat = 1;
        printf("Fail :: V_d\n");
    }

    printf("%f sec til ver2\n", (double)(clock() - t)/CLOCKS_PER_SEC);
    t = clock();
    
    //free
    free(betavals);
    free(r_tail);
    free(r_head);
    free(alpha1);
    free(alpha2);
    for(int i = 0; i < d + l; i ++)
        free(poly_b[i]);
    for(int i = 0; i < l; i ++) {
        free(poly_d[i]);
        free(poly_r[i]);
    }
    free(poly_b);
    free(poly_d);
    free(poly_r);

    //return
    return stat ? 0 : res;
}
*/

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

    mpz_init_set_str(PRIME, "1FFFFFFFFFFFFFFF", 16);
    printf("PRIME: %s\n", mpz_get_str(NULL, 10, PRIME));
  
    mpz_t* z = (mpz_t*) malloc((d + l)* sizeof(mpz_t));
    mpz_t* r = (mpz_t*) malloc((d + l) * sizeof(mpz_t));
    mpz_t* r_tail = (mpz_t*) malloc(d * sizeof(mpz_t));

    for(int i = 0; i < d + l; i ++){
//        mpz_init_set_ui(r[i], rand());
//        mod(r[i], r[i]);
        mpz_init_set_ui(r[i], 10);
    }
    for(int i = 0; i < d + l; i ++){
        mpz_init_set_ui(z[i], 10);
        //mpz_init_set_ui(z[i], rand());
        //mod(z[i], z[i]);
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
  
    evaluate_V(ri, V0, d, z);
    evaluate_V(ri0, V1, d, r_tail);
    stat = sum_check_rounding_fast(ri, C, ri, d, num_digits, num_rdigits, r, z);
    if(stat) {
        if(mpz_cmp(ri, ri0))
            printf("in the very last check, ri != ri0. ri was %s and ri0 was %s\n", mpz_get_str(NULL, 10, ri), mpz_get_str(NULL, 10, ri0));
    } else {
        printf("Failed inside the layer\n");
    }
    std::cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << std::endl;
  
   
    //set the high order of values to be those of corresponding to index i, and the low order values of z to be those corresponding to index k
    
    return 1;
}


