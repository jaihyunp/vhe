/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *Currently working: this script may not work
 *vc7.c is the latest stable script
 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*************************************
//Jai Hyun Park
//October 7, 2019.
//Implementation of fast rounding circuit using commitments
**************************************/

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define MASK 4294967295 //2^32-1
#define PRIME 2305843009213693951 //2^61-1
#define PROOF 1

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

void mod_pow(mpz_t rop, const mpz_t x, int b)
{
    int d = 1;
    mpz_t base;
    mpz_inits(base, NULL);

    mpz_set_ui(rop, 1);
    mpz_set(base, x);

    while(d <= b) {
        if(d & b)
            mod_mult(rop, rop, base);
        d = d << 1;
        if(d <= b)
            mod_mult(base, base, base);
    }

    mpz_clears(base, NULL);
}
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
    mpz_mult(rop, x, y);
    mod(rop, rop);
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
//computes b^e mod p using repeated squaring. p should be 2^61-1
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
    mpz_t v1, v2, v3, q, t1, t2, t3;
    mpz_inits(v1, v2, v3, q, t1, t2, t3, NULL);

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

    mpz_clears(v1, v2, v3, q, t1, t2, t3, NULL);
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
                mpz_set_si(tmp, r - j);
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


void mod_1neg(mpz_t rop, const mpz_t r)
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
inline void updateV(uint64* V, int num_new, uint64 ri)
{
    for(int i = 0; i < num_new; i ++) 
		V[i] = myMod(myModMult(V[i], 1 + PRIME - ri) + myModMult(V[i + num_new], ri));
}


void evaluate_V(mpz_t rop, const mpz_t *V_in, const int d, const mpz_t *r)
{
    mpz_t n, *V = (mpz_t*) malloc(sizeof(mpz_t) << d);
    mpz_init_(n);
    mpz_set_ui(n, 1);
    mpz_mult_2exp(n, n, d);
//    TODO
}
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


inline uint64 evaluate_beta(uint64 *z, uint64 *r, int d)
{
    uint64 res = 1;
    for(int i = 0; i < d; i ++)
        res = myModMult(res, myMod(myModMult(1 + PRIME - z[i], 1 + PRIME - r[i]) + myModMult(z[i], r[i])));
    return res;
}

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


/*
inline uint64 inv_pow2(int k)
{
    return myMod(1ULL << (61 - k)); 
}
*/

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

int main(int argc, char** argv)
{

    if((argc != 4)) {
        cout << "argc!= 2. Command line arg should be log(n) for nxn matrix multiplication\n";
        exit(1);
    } 
  
    int stat = 0;
    int d = atoi(argv[1]);
    int n = 1 << d;
    int num_digits = atoi(argv[2]), num_rdigits = atoi(argv[3]), del = num_digits - num_rdigits;
    int l = 0;
    while((num_digits - 1) >> ++ l);
  
    clock_t t = clock();
  
    uint64* z = (uint64*) calloc(d, sizeof(uint64));
    uint64* r = (uint64*) calloc(d + l, sizeof(uint64));
    uint64* r_tail = (uint64*) malloc(sizeof(uint64) * d);
    for(int i = 0; i < d + l; i ++)
  	    r[i] = rand() & (PRIME - 3);
//        r[i] = 3;
    for(int i = 0; i < d; i ++)
//        z[i] = 1;
  	    z[i] = rand() & (PRIME - 3);
    for(int i = 0; i < d; i ++)
        r_tail[i] = r[i + l];
    t=clock();

    //run through entire Muggles protocol with prover
    uint64* V1 = (uint64*) calloc(n, sizeof(uint64));
    uint64* V0 = (uint64*) calloc(n, sizeof(uint64));
    uint64* C = (uint64*) calloc(1 << (d + l), sizeof(uint64));
    uint64** C_old = (uint64**) malloc(num_digits * sizeof(uint64*));
    for(int i = 0; i < num_digits; i ++)
        C_old[i] = (uint64*) malloc(n * sizeof(uint64));

    t=clock();
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < num_digits; j ++){
            uint64 random_bit = rand() & 1;
            C_old[j][i] = random_bit;
            C[(i << l) + j] = random_bit;
            V1[i] ^= random_bit << j;
            if(j >= del)
                V0[i] ^= random_bit << (j - del);
        }
    }

  
    t = clock();
    uint64 ri = evaluate_V(V0, d, z);
    uint64 ri0 = evaluate_V(V1, d, r_tail);
    ri = sum_check_rounding(V1, C_old, ri, d, num_digits, num_rdigits, r_tail, z);
    if(ri != ri0) {
    	cout << "in very last check, ri != ri0. ri was : " << ri << " and ri0 was " << ri0 << endl;
    } else {
        printf("Rounding circuit verified.\n");
    }
    cout << "total OLD check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << endl;
    
  
    t = clock();
    ri = evaluate_V(V0, d, r_tail);
    ri0 = evaluate_V(V1, d, r_tail);
    ri = sum_check_rounding_fast(V1, C, ri, d, num_digits, num_rdigits, r, z);
    if(ri != ri0) {
    	cout << "in very last check, ri != ri0. ri was : " << ri << " and ri0 was " << ri0 << endl;
    } else {
        printf("Rounding circuit verified.\n");
    }
    cout << "total NEW check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << endl;
  
   
    //set the high order of values to be those of corresponding to index i, and the low order values of z to be those corresponding to index k
    
    return 1;
}


