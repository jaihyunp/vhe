/*************************************
//Justin Thaler
//February 5, 2013. 
//Implementation of Theorem 1 from the paper
//"Time-Optimal Interactive Proofs for Circuit Evaluation".
//Case Study: Matrix Multiplication
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

using namespace std;

typedef int (*fptr)(int, int);


//computes x^b
uint64 myPow(uint64 x, uint64 b)
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


//and if so subtract p, but I don't want to pay that
//efficiency hit here, so the user should just be aware of this.
uint64 myMod(uint64 x)
{
    uint64 y = x, z = y >> 61;
    while(z > 0){
        y = z + (y & PRIME);
        z = y >> 61;
    }
    if(y == PRIME)
        y = 0;
    return y;
}


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


//computes b^e mod p using repeated squaring. p should be 2^61-1
uint64 myModPow(uint64 b, uint64 e)
{
  uint64 result;
  if(e==1)
      return b;
  if(e == 0)
      return 1;
  if((e & 1) == 0) {
    result = myModPow(b, (e >> 1));
    return myModMult(result, result);
  } else {
     return myModMult(myModPow(b, e-1), b);
  }
}


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
  }while(v3 != 0 && v3!= PRIME);
}


//Computes the modular multiplicative inverse of a modulo m,
//using the extended Euclidean algorithm
//only works for p=2^61-1
uint64 inv(uint64 a)
{
  uint64 u1;
  uint64 u2;
  uint64 u3;
  extEuclideanAlg(a, &u1, &u2, &u3);

  if(u3 ==  1)
      return myMod(u1);
  else
      return 0;
}


//computes chi_v(r), where chi is the Lagrange polynomial that takes
//boolean vector v to 1 and all other boolean vectors to 0. (we view v's bits as defining
//a boolean vector. n is dimension of this vector.
//all arithmetic done mod p
uint64 chi(uint64 v, uint64* r, uint64 n)
{ 
  uint64 x=v;
  uint64 c = 1;
  for(uint64 i = 0; i <n; i++)
  {
    if(x & 1)
      c = myModMult(c, r[i]);
    else
      c = myModMult(c, 1 + PRIME - r[i]);

    x = x >> 1;
  }
  return c;
}

//evaluates beta_z polynomial (described in GKR08)
//at location k (the bits of k are interpreted as a boolean vector). mi is dimension of k and r.
uint64 evaluate_beta_z_ULL(int mi, uint64* z, uint64 k)
{
  uint64 ans = 1;
  uint64 x = k;
  for(int i = 0; i < mi; i ++)
  {
    if(x & 1)
      ans = myModMult(ans, z[i]);
    else
      ans = myModMult(ans, 1 + PRIME - z[i]);
    x = x >> 1;
  }
  return ans;
}


//evaluates V_i polynomial at location r.
//Here V_i is described in GKR08; it is the multi-linear extension
//of the vector of gate values at level i of the circuit
uint64 evaluate_V_i(int mi, int ni, uint64* level_i, uint64* r)
{
  uint64 ans = 0;
  for(uint64 k = 0; k < ni; k ++)
  {
     ans=myMod(ans + myModMult(level_i[k], chi(k, r, mi)));
  }
  return ans;
}


//extrapolate the polynomial implied by vector vec of length n to location r
uint64 extrap(uint64* vec, uint64 n, uint64 r)
{
  uint64 result = 0;
  uint64 mult = 1;
  for(uint64 i = 0; i < n; i ++)
  {
    mult = 1;
    for(uint64 j = 0; j < n; j ++)
    {
      if (i > j)
        mult = myModMult(myModMult(mult, myMod(r - j + PRIME)), inv(i - j));
      if (i < j)
        mult = myModMult(myModMult(mult, myMod(r - j + PRIME)), inv(myMod(i + PRIME - j)));
    }
    result = myMod(result + myModMult(mult, vec[i]));
  }
  return result;
}

void initialize_alphavals(uint64* alphavals, int min, int max);
void update_alphavals(uint64* alphavals, int num_new, uint64 rval);
void update_V_b();
uint64 sum_check_rounding_fast(~)
{
}

void updateV(uint64* V, int num_new, uint64 ri)
{
    for(int i = 0; i < num_new; i ++) {
		V[i] = myMod(myModMult(V[i], 1 + PRIME - ri) + myModMult(V[i + num_new], ri));
	}
}


//sets betavals(p) = \prod_{i=0}^{d-1} (zipi + (1-zi)(1-pi))
void initialize_betavals(uint64* betavals, int d, uint64* z)
{
	uint64 zval;
	int two_to_k = 1;
	uint64 oldval;
	
    betavals[0] = 1;
	for(int k = 0; k < d; k++)
	{
		zval = z[k];
		for(int i = 0; i < two_to_k; i++)
		{
			oldval = betavals[i];
			betavals[i] = myModMult(oldval, 1 + PRIME - zval);
			betavals[i + two_to_k] = myModMult(oldval, zval);
		}
		two_to_k = two_to_k * 2;
	}
}


void update_betavals(uint64* betavals, int num_new, uint64 rval)//이렇게 하면 굳이 zval을 받을 필요 없고, inv할 필요없음.
{
    for(int i = 0; i < num_new; i ++){
        betavals[i] = myMod(myModMult(betavals[i], 1 + PRIME - rval) + myModMult(betavals[i + num_new], rval));
    }
}


uint64 sum_check_rounding(uint64* V, uint64 **C, uint64 ri, int d, int num_digits, int num_rdigits, uint64 *r, uint64 *z)
{
    uint64 ***poly, *betavals;
    uint64 stat = 0;

    //Assumption: Commitment enables us to evaluate C[i] with verification
    uint64 *cr = (uint64*) malloc(num_digits * sizeof(uint64));
    uint64 *cz = (uint64*) malloc(num_digits * sizeof(uint64));
    for(int i = 0; i < num_digits; i ++){
        cr[i] = evaluate_V_i(d, 1 << d, C[i], r);
        cz[i] = evaluate_V_i(d, 1 << d, C[i], z);
    }


    //initialize poly
    poly = (uint64***) malloc(num_digits * sizeof(uint64**));
    for(int i = 0; i < num_digits; i ++){
        poly[i] = (uint64**) malloc(d * sizeof(uint64*));
        for(int j = 0; j < d; j ++){
            poly[i][j] = (uint64*) calloc(4, sizeof(uint64));
        }
    }

    //initialize betavals
    betavals = (uint64*) malloc((1 << d) * sizeof(uint64));
    initialize_betavals(betavals, d, z);


    //prover do:
    int num_terms =  1 << d;
    for(int round = 0; round < d; round ++){

        uint64 tmp = num_terms >> 1;
        for(int i = 0; i < num_digits; i ++){
            for(int j = 0; j < num_terms / 2; j ++){

                uint64 b0 = betavals[j],
                       b1 = betavals[j + tmp],
                       b2 = myMod(2 * b1 + PRIME - b0),
                       b3 = myMod(3 * b1 + 2 * 2 * PRIME - 2 * b0);
                uint64 c0 = C[i][j],
                       c1 = C[i][j + tmp],
                       c2 = myMod(2 * c1 + PRIME - c0),
                       c3 = myMod(3 * c1 + 2 * PRIME - 2 * c0);
                
                poly[i][round][0] = myMod(poly[i][round][0] + myModMult(myModMult(b0, c0), c0 + PRIME - 1));
                poly[i][round][1] = myMod(poly[i][round][1] + myModMult(myModMult(b1, c1), c1 + PRIME - 1));
                poly[i][round][2] = myMod(poly[i][round][2] + myModMult(myModMult(b2, c2), c2 + PRIME - 1));
                poly[i][round][3] = myMod(poly[i][round][3] + myModMult(myModMult(b3, c3), c3 + PRIME - 1));
            }

    		updateV(C[i], num_terms >> 1, r[d - round - 1]);
        }
            
        num_terms >>= 1;
		update_betavals(betavals, num_terms, r[d - round - 1]);
    }


    //verifier do:
    uint64 vr = 0;
    for(int i = num_digits - num_rdigits; i < num_digits; i ++)
        vr = myMod(vr + myModMult(1ULL << (i - num_digits + num_rdigits), cz[i]));
    if(vr != ri){
        printf("Fail :: not match with the previous layer\n");
        stat = 1;
    } else {
        printf("Success :: prev layer\n");
    }
    
    if(!stat) {
        for(int i = 0; i < num_digits; i ++){
            if(myMod(poly[i][0][0] + poly[i][0][1]) != 0){
                printf("Fail :: %dth digit is_bit sumcheck init: %lld should be zero\n", i, poly[i][0][1] + poly[i][0][1]);
                stat = 1;
                break;
            } else { 
                printf("Success :: %dth digit is_bit 1st sumcheck\n", i);
            }

            uint64 extrap_val = extrap(poly[i][0], 4, r[d - 1]);
            for(int round = 1; round < d; round ++){
                if(myMod(poly[i][round][0] + poly[i][round][1]) != extrap_val){
                    printf("Fail :: %dth digit is_bit sumcheck %dth round: %lld should be %lld\n", i, round, myMod(poly[i][round][0] + poly[i][round][1]), extrap_val);
                    stat = 1;
                    break;
                } else {
                    printf("Success :: %dth digit is_bit sumcheck %dth round\n", i, round);
                    extrap_val = extrap(poly[i][round], 4, r[d - round - 1]);
                }
            }

            if(!stat){
                if(extrap_val != myModMult(betavals[0], myModMult(cr[i], cr[i] + PRIME - 1))){
                    printf("Fail :: %dth digit is_bit sumcheck does not match with the commitment\n", i);
                    stat = 1;
                } else {
                    printf("Success :: %dth digit is_bit sumcheck commitment\n", i);
                }
            }
        }
    }
    if(!stat){
        for(int i = 0; i < num_digits; i ++)
            stat = myMod(stat + myModMult(1ULL << i, C[i][0]));
    }


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
    return stat;
}
/*
//v1 contains MLE under x0 = x0_(d-1) x0_d ... x0_0, x1, ..., xn
//v2 contains MLE under xi_j * (2 ^ j)
uint64 sum_check_rounding(uint64* V, uint64 *V0, uint64 *V1, uint64 ri, int log_num_elems, int num_digits, int num_rdigits, uint64* r, uint64 *z)
{
    int stat = 0;
	int num_elems = 1 << log_num_elems;
    int up_num_digits = 1, log_num_digits = 0;
    while(up_num_digits < num_digits){
        up_num_digits <<= 1;
        log_num_digits ++;
    }
    int log_num_bins = log_num_elems + log_num_digits, num_bins = 1 << log_num_bins;


    uint64 **poly_b0, **poly_b1, **poly_d, **poly_r; //poly_b for is_bit, poly_d for is_digit, poly_r for rounding
    uint64 *betavals_b, *betavals_d;

    uint64 *v3 = (uint64*) malloc(up_num_digits * sizeof(uint64));
    uint64 *v4 = (uint64*) malloc(up_num_digits * sizeof(uint64));

    //initialize v3 & v4
    for(int i = 0; i < up_num_digits; i ++){
        if(i < num_rdigits){
            v3[i] = 1 << i;
            v4[i] = 1 << i;
        } else if (i < num_digits) {
            v3[i] = 1 << i;
            v4[i] = 0;
        } else {
            v3[i] = 0;
            v4[i] = 0;
        }
    }

    //initialize poly
    poly_b0 = (uint64**) malloc(log_num_bins * sizeof(uint64*));
    poly_b1 = (uint64**) malloc(log_num_bins * sizeof(uint64*));
    poly_d = (uint64**) malloc(log_num_bins * sizeof(uint64*));
    poly_r = (uint64**) malloc(log_num_bins * sizeof(uint64*));
    for(int i = 0; i < log_num_elems; i ++){
        poly_d[i] = (uint64*) calloc(4, sizeof(uint64));
        poly_r[i] = (uint64*) calloc(4, sizeof(uint64));
    }
    for(int i = 0; i < log_num_bins; i ++){
        poly_b0[i] = (uint64*) calloc(4, sizeof(uint64));
        poly_b1[i] = (uint64*) calloc(4, sizeof(uint64));
    }

    
    //initialize betavals
    betavals_b = (uint64*) malloc(num_bins * sizeof(uint64));
    betavals_d = (uint64*) malloc(num_elems * sizeof(uint64));
	initialize_betavals(betavals_b, log_num_elems + log_num_digits, z_b);
	initialize_betavals(betavals_d, log_num_elems, z_d);


    int num_terms = num_elems;
	for(int round = 0; round < log_num_elem; round ++) {

        for(int j = 0; j < num_terms / 2; j ++) {

            uint64 tmp = num_terms >> 1;

            //betavals for poly_b
            uint64 bb0 = betavals_b[j],
                   bb1 = betavals_b[j + tmp],
                   bb2 = myMod(2 * bb1 + PRIME - bb0),
                   bb3 = myMod(3 * bb1 + PRIME - 2 * bb0);

            //v vals for poly_b
            uint64 vb0 = myModMult(V2[j], myMod(V2[j] + PRIME - 1)),
                   vb1 = myModMult(V2[j + tmp], myMod(V2[j + tmp] + PRIME - 1)),
                   vb2 = myMod(2 * vb1 + PRIME - vb0),
                   vb3 = myMod(3 * vb1 + PRIME - vb0);

            //v vals for poly_d
            uint64 vd0 = myMod(PRIME - V1[j]),
                   vd1 = myMod(PRIME - V1[j + tmp]);
            for(int i = 0; i < num_digits; i ++) {
                vd0 = myMod(vd0 + myModMult(1 << i, V2[(j << num_digits) + i]));
                vd1 = myMod(vd1 + myModMult(1 << i, V2[(j << num_digits) + i + tmp]));
            }
            uint64 vd2 = myMod(2 * vd1 + PRIME - vd0);

            //v vals for poly_r
//            uint64 vr0 = 0,
//                   vr1 = 0;
//            for(int i = 0; i < num_rdigits; i ++) {
//                vr0 = myMod(vr0 + myModMult(1 << i, V2[(j << num_digits) + i]));
//                vr1 = myMod(vr1 + myModMult(1 << i, V2[(j << num_digits) + i + tmp]));
//            }
//            uint64 vr2 = myMod(2 * vr1 + PRIME - vr0);


            //evaluate polys at 0, 1 and 2
            poly_b[round][0] = myMod(poly_b[round][0] + myModMult(bb0, vb0));
            poly_b[round][1] = myMod(poly_b[round][1] + myModMult(bb1, vb1));
            poly_b[round][2] = myMod(poly_b[round][2] + myModMult(bb2, vb2));

            poly_d[round][0] = myMod(poly_d[round][0] + myModMult(bd0, vd0));
            poly_d[round][1] = myMod(poly_d[round][1] + myModMult(bd1, vd1));
            poly_d[round][2] = myMod(poly_d[round][2] + myModMult(bd2, vd2));

//            poly_r[round][0] = myMod(poly_r[round][0] + myModMult(br0, vr0));
//            poly_r[round][1] = myMod(poly_r[round][1] + myModMult(br1, vr1));
//            poly_r[round][2] = myMod(poly_r[round][2] + myModMult(br2, vr2));
        }

		num_terms = num_terms >> 1;
        
		update_betavals(betavals_b, num_terms, r[log_bins - round]);
		update_betavals(betavals_d, num_terms, r[log_bins - round]);

		updateV(V1, num_terms, r[log_bins - round]);
        updateV(V2, num_terms << (log_num_digits), r[log_bins - round]);
        update
	}

	

    //From here is what Verifier do:
    if(myMod(poly_b[0][0] + poly_b[0][1]) != 0) {
        printf("not bit\n");
        stat = 1;
    } else if(myMod(poly_d[0][0] + poly_d[0][1]) != 0) {
        printf("not digit\n");
        stat = 1;
    } else if(myMod(poly_r[0][0] + poly_r[0][1]) != ri) {
        printf("1st check sum rounding\n");
        stat = 1;
    } else {
        uint64 extrap_val_b = extrap(poly_b[0], 3, r[log_bins]);
        uint64 extrap_val_d = extrap(poly_d[0], 3, r[log_bins]);
        uint64 extrap_val_r = extrap(poly_r[0], 3, r[log_bins]);

        for(int round = 1; roumd < log_elems; round ++){
            if(myMod(poly_b[round][0] + poly_b[round][1]) != extrap_val_b) {
                printf("check for bit %d round fails\n", round);
                stat = 1;
                break;
            } else if(myMod(poly_d[round][0] + poly_d[round][1]) != extrap_val_d) {
                printf("check for digit %d round fails\n", round);
                stat = 1;
                break;
            } else if(myMod(poly_r[round][0] + poly_r[round][1]) != extrap_val_r) {
                pritnf("check for round %d round fails\n", round);
                stat = 1;
                break;
            } else {
                uint64 extrap_val_b = extrap(poly_b[0], 3, r[log_bins]);
                uint64 extrap_val_d = extrap(poly_d[0], 3, r[log_bins]);
                uint64 extrap_val_r = extrap(poly_r[0], 3, r[log_bins]);
            }
        }
    }

	if(myMod(poly[0][0] + poly[0][1]) != 0)	{
		cout << "first check failed when ni was " << ni << endl;
    	stat = 1;
	} else {
        uint64 extrap_val = extrap(poly[0], 3, r[log_elems]);
        for(int round = 1; round < log_elems; round++)
        {
            if(myMod(poly[round][0] + poly[round][1]) != extrap_val) 
            {
                cout << "check for round " << round << " failed when ni was " << ni << endl;
                cout << "extrap_val was " << extrap_val << endl;
                cout << "poly[round][0] was " << poly[round][0] << " and poly[round][1] was " << poly[round][1] << " and their sum was " << myMod(poly[round][0] + poly[round][1]) << endl;
                stat = 1;
                break;
            } else {
                extrap_val = extrap(poly[round], 3, r[d - round]);
            }
         }
    }         
	 //now p "tells" v \tilde{v}(p)  and v checks final message based on this
	 //i.e. v checks g_{d-1}(r) = beta(z, p) \tilde{v}(rd, ..., r1, r0). 
	 //if so, v believes p as long as \tilde{v}(p) = extrap(v, 2, r)
	
    if(!stat){ 
    	uint64 correct_out = myModMult(betavals[0], myMod(v[0] + v[1]));
	    if(correct_out != extrap_val) {
	    	cout << "correct_out != extrap_val. correct_out is " << correct_out << " and extrap_val is: " << extrap_val << endl;
            stat = 1;
    	}
    }

    if(!stat)
        stat = extrap(V, 2, r[d]);
    else
        stat = -1;

	return stat;
}
*/
//NO USE
/*
uint64 sum_check_rounding(uint64* v1, uint64 *v2, uint64 ri, int log_num_elems, int log_num_digits, uint64* r1, uint64 *r2, uint64* z, uint64 *z_b, uint64 *z_d)
{
	int num_elems = 1 << log_num_elems;
    int num_digits = 1 << log_num_digits;
    int num_bins = 1 << (log_num_elems + log_num_digits);
    int stat = 0;

    uint64 **poly_b, **poly_d, **poly_r; //poly_b for is_bit, poly_d for is_digit, poly_r for rounding
    uint64 *betavals;

    
    //initialize poly
    poly_b = (uint64**) malloc((log_num_elems + log_num_digits) * sizeof(uint64*));
    poly_d = (uint64**) malloc(log_num_elems * sizeof(uint64*));
    poly_r = (uint64**) malloc(log_num_elems * sizeof(uint64*));
    for(int i = 0; i < log_num_elems; i ++) {
        poly_d[i] = (uint64*) calloc(4, sizeof(uint64));
        poly_r[i] = (uint64*) calloc(4, sizeof(uint64));
    }
    for(int i = 0; i < log_num_elems + log_num_digits; i ++)
        poly_b[i] = (uint64*) calloc(4, sizeof(uint64));

    //initialize betavals
	initialize_betavals(betavals, log_num_elems, z);//d-1 degree poly * 2


	for(int round = 0; round < d; round++) {
        for(int j = 0; j < num_terms / 2; j ++) {
            uint64 tmp = num_terms >> 1;
            uint64 b0 = betavals[j],
                   b1 = betavals[j + tmp],
                   b2 = myMod(2 * b1 + PRIME - b0);
            uint64 v0 = myModMult(v[j], myMod(v[j] + PRIME - 1)),
                   v1 = myMod(v[j + tmp], myMod(v[j + tmp] + PRIME - 1),
                   v2 = myMod(2 * v1 + PRIME - v0);

            poly[round][0] = myMod(poly[round][0] + myModMult(b0, v0));
            poly[round][1] = myMod(poly[round][1] + myModMult(b1, v1));
            poly[round][2] = myMod(poly[round][2] + myModMult(b2, v2));
        }

		num_terms = num_terms >> 1;
        
		update_betavals(betavals, num_terms, 0, r[d - round]);
		updatev(v, num_terms << 1, r[d - round]);
	}
	
	if(myMod(poly[0][0] + poly[0][1]) != 0)	{
		cout << "first check failed when ni was " << ni << endl;
    	stat = 1;
	} else {
        uint64 extrap_val = extrap(poly[0], 3, r[d]);
        for(int round = 1; round < d; round++)
        {
            if(myMod(poly[round][0] + poly[round][1]) != extrap_val) 
            {
                cout << "check for round " << round << " failed when ni was " << ni << endl;
                cout << "extrap_val was " << extrap_val << endl;
                cout << "poly[round][0] was " << poly[round][0] << " and poly[round][1] was " << poly[round][1] << " and their sum was " << myMod(poly[round][0] + poly[round][1]) << endl;
                stat = 1;
                break;
            } else {
                extrap_val = extrap(poly[round], 3, r[d - round]);
            }
         }
    }         
	 //now p "tells" v \tilde{v}(p)  and v checks final message based on this
	 //i.e. v checks g_{d-1}(r) = beta(z, p) \tilde{v}(rd, ..., r1, r0). 
	 //if so, v believes p as long as \tilde{v}(p) = extrap(v, 2, r)
	
    if(!stat){ 
    	uint64 correct_out = myModMult(betavals[0], myMod(v[0] + v[1]));
	    if(correct_out != extrap_val) {
	    	cout << "correct_out != extrap_val. correct_out is " << correct_out << " and extrap_val is: " << extrap_val << endl;
            stat = 1;
    	}
    }

    if(!stat)
        stat = extrap(v, 2, r[d]);
    else
        stat = -1;

	return stat;
}

*/




////////////////////////////////TODO/////////////////////////////////
/*
void updateV2(uint64* V, int num_new, int log_num_digits, uint64 *r2, int ind)
{
    int num_terms = num_new * (1 << log_num_digits), num_terms0 = (ind + 1) * log_num_digits;

    for(int j = 0; j < log_num_digits; j ++) {

        for(int i = 0; i < num_terms; i ++)
            V[i] = myMod(myModMult(V[i], 1 + PRIME - r2[num_terms0 - j]) + myModMult(V[i + num_terms], r2[num_terms0 - j]));

        num_terms >>= 1;
    }
}


uint64 sum_check_is_binary(uint64* V1, uint64* V2, uint64 ri, int log_num_digits, int num_elem, uint64* r1, uint64* r2, uint64** poly, uint64* betavals, uint64* z)
{

	for(int i = 0; i < d; i ++)
		poly[i][0] = poly[i][1] = poly[i][2] = poly[i][3] = 0;
	
	initialize_betavals(betavals, d, z);//d-1 degree poly * 2
	int num_terms = ni;//n

	for(int round = 0; round < d; round ++)
	{
        for(int j = 0; j < num_terms / 2; j ++){
            uint64 tmp = num_terms >> 1;
            uint64 b0 = betavals[j],
                   b1 = betavals[j + tmp],
                   b2 = myMod(2 * b1 + PRIME - b0);
            uint64 v0 = myModMult(V[j], myMod(V[j] + PRIME - 1)),
                   v1 = myMod(V[j + tmp], myMod(V[j + tmp] + PRIME - 1),
                   v2 = myMod(2 * v1 + PRIME - v0);

            poly[round][0] = myMod(poly[round][0] + myModMult(b0, v0));
            poly[round][1] = myMod(poly[round][1] + myModMult(b1, v1));
            poly[round][2] = myMod(poly[round][2] + myModMult(b2, v2));

        }

		num_terms = num_terms >> 1;
        
		update_betavals(betavals, num_terms, 0, r1[d - round]);
		updateV(V, num_terms << 1, r1[d - round]);
        updateV2(V, num_terms, log_num_digis, r2[d - round], round); 
	}
	
	if(myMod(poly[0][0] + poly[0][1]) != 0)
	{
		cout << "first check failed when ni was " << ni << endl;
		cout << "ri was " << 0 << endl;
		cout << "poly[0][0] was " << poly[0][0] << " and poly[0][1] was " << poly[0][1] << " and their sum was " << myMod(poly[0][0] + poly[0][1]) << endl;
		exit(1);
	}

	uint64 extrap_val = extrap(poly[0], 3, r[d]);
	for(int round = 1; round < d; round++)
	{
		if(myMod(poly[round][0] + poly[round][1]) != extrap_val) 
		{
			cout << "check for round " << round << " failed when ni was " << ni << endl;
			cout << "extrap_val was " << extrap_val << endl;
			cout << "poly[round][0] was " << poly[round][0] << " and poly[round][1] was " << poly[round][1] << " and their sum was " << myMod(poly[round][0] + poly[round][1]) << endl;
			exit(1);
	 	}
	 	extrap_val = extrap(poly[round], 3, r[d - round]);
	 }
	 
	 //now P "tells" V \tilde{V}(p)  and V checks final message based on this
	 //i.e. V checks g_{d-1}(r) = beta(z, p) \tilde{V}(rd, ..., r1, r0). 
	 //If so, V believes P as long as \tilde{V}(p) = extrap(V, 2, r)
	 
	 uint64 correct_out = myModMult(betavals[0], myMod(V[0] + V[1]));
	 if(correct_out != extrap_val)
	 {
	 	cout << "correct_out != extrap_val. correct_out is " << correct_out << " and extrap_val is: " << extrap_val << endl;
	 	exit(1);
	 }
	 return extrap(V, 2, r[d]);

}
*/
uint64 sum_check_polyadd(uint64* V, uint64 ri, int d, int ni, int* com_ct, int* rd_ct, uint64* r, uint64** poly, uint64* betavals, uint64* z)
{

	for(int i = 0; i < d + 1; i ++){
		poly[i][0] = poly[i][1] = poly[i][2] = poly[i][3] = 0;
    }

	*com_ct = *com_ct + 3 * d + 2;
	*rd_ct = *rd_ct+d + 1;
	
	initialize_betavals(betavals, d, z);//d-1 degree poly * 2
	int num_terms = ni;//n

	for(int round = 0; round < d; round++)
	{
        for(int j = 0; j < num_terms / 2; j ++){
            uint64 tmp = num_terms >> 1;
            uint64 b0 = betavals[j],
                   b1 = betavals[j + tmp],
                   b2 = myMod(2 * b1 + PRIME - b0);
            uint64 v0 = myMod(V[j] + V[j + num_terms]),
                   v1 = myMod(V[j + tmp] + V[j + tmp + num_terms]),
                   v2 = myMod(2 * v1 + PRIME - v0);

            poly[round][0] = myMod(poly[round][0] + myModMult(b0, v0));
            poly[round][1] = myMod(poly[round][1] + myModMult(b1, v1));
            poly[round][2] = myMod(poly[round][2] + myModMult(b2, v2));

        }

		num_terms = num_terms >> 1;
        
		update_betavals(betavals, num_terms, r[d - round - 1]);
		updateV(V, num_terms << 1, r[d - round - 1]);
	}
	
	if( (myMod(poly[0][0] + poly[0][1]) != ri) && (myMod(poly[0][0] + poly[0][1]) != ri + PRIME) && 
		((myMod(poly[0][0] + poly[0][1]) + PRIME) != ri))
	{
		cout << "first check failed when ni was " << ni << endl;
		cout << "ri was " << ri << endl;
		cout << "poly[0][0] was " << poly[0][0] << " and poly[0][1] was " << poly[0][1] << " and their sum was " << myMod(poly[0][0] + poly[0][1]) << endl;
		exit(1);
	}

	uint64 extrap_val = extrap(poly[0], 3, r[d - 1]);
	for(int round = 1; round < d; round++)
	{
		if( (myMod(poly[round][0] + poly[round][1]) != extrap_val) && (myMod(poly[round][0] + poly[round][1]) != extrap_val + PRIME) &&
			((myMod(poly[round][0] + poly[round][1]) + PRIME) != extrap_val)) 
		{
			cout << "check for round " << round << " failed when ni was " << ni << endl;
			cout << "extrap_val was " << extrap_val << endl;
			cout << "poly[round][0] was " << poly[round][0] << " and poly[round][1] was " << poly[round][1] << " and their sum was " << myMod(poly[round][0] + poly[round][1]) << endl;
			exit(1);
	 	}
	 	extrap_val = extrap(poly[round], 3, r[d - round - 1]);
	 }
	 
	 //now P "tells" V \tilde{V}(p)  and V checks final message based on this
	 //i.e. V checks g_{d-1}(r) = beta(z, p) \tilde{V}(rd, ..., r1, r0). 
	 //If so, V believes P as long as \tilde{V}(p) = extrap(V, 2, r)
	 
	 uint64 correct_out = myModMult(betavals[0], myMod(V[0] + V[1]));
	 if( (correct_out != extrap_val) && (correct_out + PRIME != extrap_val) && (correct_out != extrap_val + PRIME))
	 {
	 	cout << "correct_out != extrap_val. correct_out is " << correct_out << " and extrap_val is: " << extrap_val << endl;
	 	exit(1);
	 }
	 return extrap(V, 2, r[d]);
}
	

int main(int argc, char** argv)
{

  if((argc != 4)) {
    cout << "argc!= 2. Command line arg should be log(n) for nxn matrix multiplication\n";
    exit(1);
  }
  
  int d = atoi(argv[1]);
  int n = 1 << d;
  clock_t t = clock();
  
  uint64* z = (uint64*) calloc(d, sizeof(uint64));
  uint64* r = (uint64*) calloc(d, sizeof(uint64));

  for(int i = 0; i < d + 1; i ++)
  	r[i] = rand() + 4;
  for(int i = 0; i < d; i ++)
    z[i] = rand() + 4;

  t=clock();

  //run through entire Muggles protocol with prover
  int num_digits = atoi(argv[2]), num_rdigits = atoi(argv[3]), del = num_digits - num_rdigits;
  uint64* V1 = (uint64*) calloc(n, sizeof(uint64));
  uint64* V0 = (uint64*) calloc(n, sizeof(uint64));
  uint64** C = (uint64**) malloc(num_digits * sizeof(uint64*));
  for(int i = 0; i < num_digits; i ++)
      C[i] = (uint64*) malloc(n * sizeof(uint64));

  t=clock();
  for(int i = 0; i < num_digits; i ++){
        for(int j = 0; j < n; j ++){
            C[i][j] = rand() & 1;
            V1[j] ^= C[i][j] << i;
            if(i >= del)
                V0[j] ^= C[i][j] << (i - del);
        }
  }

    printf("Initialization:\n");
    int check;
    for(int i = 0; i < n; i ++){
        printf("%3d:", i);
        check = 0;
        for(int j = 60; j >= 0; j -= 4) {
            if(check || ((V1[i] >> j) & 0xF)){
                printf("%X", (unsigned int)((V1[i] >> j) & 0xF));
                check = 1;
            }
        }
        printf("\n -> ");
        check = 0;
        for(int j = 60; j >= 0; j -= 4){
            if(check || ((V0[i] >> j) & 0xF)){
                printf("%X", (unsigned int)((V0[i] >> j) & 0xF));
                check = 1;
            }
        }
        printf("\n");
    }

 t = clock();
  uint64 ri = evaluate_V_i(d, n, V0, z);
  uint64 ri0 = evaluate_V_i(d, n, V1, r);
//  ri = sum_check_polyadd(V1, ri, d, n, &com_ct, &rd_ct, r, poly, betavals, z);
  ri = sum_check_rounding(V1, C, ri, d, num_digits, num_rdigits, r, z);


  cout << "total check time is: " << (double)((double) clock()-t)/CLOCKS_PER_SEC << endl;
 
  //set the high order of values to be those of corresponding to index i, and the low order values of z to be those corresponding to index k
  
  t=clock();	
  if(ri != ri0) {
  	cout << "in very last check, ri != ri0. ri was : " << ri << " and ri0 was " << ri0 << endl;
  	//exit(1);
  } else {
      printf("Rounding circuit verified.\n");
  }

  return 1;
}

