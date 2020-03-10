<!---Implementation of GKR Protocol on (somewhat) HE-->
# Implementation of VC over Cyclotomic Rings.

## Description of the structure of each experiments:

### Multiplication over F[X].

Suppose that we are given **num** pairs of **N** degree polynomials, f~i~ and g~i~ in F[X], where F is a finite field of order prime. Vf delegates the evaluation of a circuit computing h~i~ = f~i~ g~i~.
The entire process would be as follows.

1. V3\[\]\[\] is a double array of F, in which V3\[2\*i\] contains the coefficients of f~i~ and V\[2\*i\] contains that of g~i~. Vf sends her V3 to Pv. 

2. Pv computes a double array of F V2\[\]\[\], which contains the coefficients of f~i~ g~i~. Pv then sends V2 to Vf.

3. Vf picks a random point **val** from F, and evaluates the polynomials in input and output layers at val. She sends val to Pv, and Pv evaluates all polynomials in the circuit at val; this would be represented as a *multi-linear map*.

4. The GKR Protocol. Exploit the equation:
$$
V2(z) = \sum \beta(z;p) V3(p0) V3(p1).
$$
Pv may send the GKR polynomials to Vf with auxiliary values: V3(r0) and V3(r1). Vf can check whether the output of GKR protocol is consistent to V3(r\*). Then she can reduce the output to (expected) V3(r). 

5. Vf checks whether the expected value of V3 is consistent to the original V3.



### Multiplication over F[X]/(X^N +1).

Suppose that we are given **num** pairs of elements in F[x]/(X^N +1), f~i~ and g~i~, where F is a finite field of order prime. Vf delegates the evaluation of a circuit computing h~i~ = f~i~ g~i~.
The entire process would be as follows.

1. V3\[\]\[\] is a double array of F, in which V3\[2\*i\] contains the coefficients of the representative element of f~i~ and V\[2\*i\] contains that of g~i~. Vf sends her V3 to Pv. 
**Vf additionally sends the Commit keys to Pv**.

2. Pv computes a double array of F V2\[\]\[\], which contains the coefficients of f~i~ g~i~, **and commit the coefficients**; denote it C1. He reduces the result as each polynomial has degree less than N; and denote it as V1. Pv sends V1 to Vf.

3. Vf picks a random point **val** from F, and evaluates the polynomials in input (V3) and output (V1) layers at val. She sends val to Pv, and Pv evaluates all polynomials in the circuit (V3, V2, and V1) at val; this would be represented as a *multi-linear map*.

4. The GKR Protocol. Exploit the equation:
$$
(poly) V2(z) = \sum \beta(z;p) V3(p0) V3(p1).

(polyC10) V2(z) = \sum \beta(z;p) \tau(q) (-C1(p,q,1)+C1(p,q,0)).

(polyC11) V1(z) = \sum \beta(z;p) \tau(q) (C1(p,q,1)val^N^ +C1(p,q,0)).
$$
Pv may send the GKR polynomials to Vf with auxiliary values: V3(r0), V3(r1), C1(r,0), and C1(r,1). Vf can check whether the output of GKR protocol is consistent to V3(r\*) and C1(r\*). Then she can reduce the output to (expected) V3(r) and C1(r). 

5. Vf checks whether the expected values of V3 and C1(r) is consistent to the original V3 and committed C1.




---
## Choice of PRIME

PRIME is the prime number, which denotes the field that all VC computations play. In order to compute efficiently, we'd better nicely choose PRIME as follows:
1. In order to find the (4N)th root of unity efficiently, I pick PRIME as a prime number of form 4Nk+1, in which k is a prime number as well.
2. For the use of commitment, PRIME^2 + 1 (= GSIZE) should be a prime number as well. (NOT IMPLEMENTED YET)


## The Use of ROU

ROU stores the cyclic group generated by a 4N-th root of unity.
1. This is for FFT 

