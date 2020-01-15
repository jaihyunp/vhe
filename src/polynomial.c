#include "parameters.h"
#include "field.h"
#include "polynomial.h"
#include <stdlib.h>
#include <gmp.h>

void _fft(mpz_t *buf, mpz_t *out, const uint64 n, const uint64 step)
{
    mpz_t tmp1; 
    mpz_init(tmp1);

    if(step < n) {
        _fft(out, buf, n, step * 2);
        _fft(out + step, buf + step, n, step * 2);

        for(uint64 i = 0; i < n; i += 2 * step) {
            mod_mult(tmp1, ROU[i * (2 * N / n)], out[i + step]);//2n-th rou * i
            mod_add(buf[i / 2], out[i], tmp1);
            mod_sub(buf[(i + n) / 2], out[i], tmp1);
        }
    }

    mpz_clear(tmp1);
}

void fft(mpz_t *v, const uint64 n, const mpz_t *c, const uint64 deg)
{
    mpz_t *v1 = (mpz_t*) malloc(sizeof(mpz_t) * n);
    
    for(uint64 i = 0; i < n; i ++) {
        if(i < deg){
            mpz_init_set(v1[i], c[i]);
            mpz_set(v[i], c[i]);
        } else {
            mpz_init_set_ui(v1[i], 0);
            mpz_set_ui(v[i], 0);
        }
    }
    _fft(v, v1, n, 1);

    for(uint64 i = 0; i < n; i ++)
        mpz_clear(v1[i]);
    free(v1);
}

void _ifft(mpz_t *buf, mpz_t *out, const uint64 n, const uint64 step)
{
    mpz_t tmp1; 
    mpz_init(tmp1);

    if(step < n) {
        _ifft(out, buf, n, step * 2);
        _ifft(out + step, buf + step, n, step * 2);

        for(uint64 i = 0; i < n; i += 2 * step) {
            if(i == 0)
                mod_mult(tmp1, ROU[i], out[i + step]);
            else
                mod_mult(tmp1, ROU[4 * N - i * (2 * N / n)], out[i + step]);
            mod_add(buf[i / 2], out[i], tmp1);
            mod_sub(buf[(i + n) / 2], out[i], tmp1);
        }
    }

    mpz_clear(tmp1);
}

void ifft(mpz_t *c, const mpz_t *v, const uint64 n)
{
    mpz_t *c1 = (mpz_t*) malloc(sizeof(mpz_t) * n), tmp;
    mpz_init(tmp);

    for(uint64 i = 0; i < n; i ++) {
        mpz_init_set(c1[i], v[i]);
        mpz_set(c[i], v[i]);
    }
    _ifft(c, c1, n, 1);

    mpz_set_ui(tmp, n);
    mod_inv(tmp, tmp);
    for(uint64 i = 0; i < n; i ++) {
        mod_mult(c[i], c[i], tmp);//c[i]<-c[i]/N
    }

    mpz_clear(tmp);
    for(uint64 i = 0; i < n; i ++)
        mpz_clear(c1[i]);
    free(c1);
}


void coeff_mult(mpz_t *c, const mpz_t *c1, const mpz_t *c2, const uint64 deg)
{
    mpz_t tmp;
    mpz_init(tmp);
    for(uint64 i = 0; i < 2 * deg; i ++) {
        mpz_set_ui(c[i], 0);
    }
    for(uint64 i = 0; i < deg; i ++) {
        for(uint64 j = 0; j < deg; j ++) {
            mod_mult(tmp, c1[i], c2[j]);
            mod_add(c[i + j], c[i + j], tmp);
        }
    }
    mpz_clear(tmp);
}

void fourier_mult(mpz_t *v, const mpz_t *v1, const mpz_t *v2, const uint64 n)
{
    for(uint64 i = 0; i < n; i ++)
        mod_mult(v[i], v1[i], v2[i]);
}

void coeff_add(mpz_t *c, const mpz_t *c1, const mpz_t *c2, const uint64 deg)
{
    for(uint64 i = 0; i < deg; i ++)
        mod_add(c[i], c1[i], c2[i]);
}

void fourier_add(mpz_t *v, const mpz_t *v1, const mpz_t *v2, const uint64 n)
{
    for(uint64 i = 0; i < n; i ++)
        mod_add(v[i], v1[i], v2[i]);
}

int polynomial_cmp(const mpz_t *c, const uint64 deg, const mpz_t *v, const uint64 n)
{
    int stat = 1;
    mpz_t *tmp = (mpz_t*) malloc(sizeof(mpz_t) * n);
    for(uint64 i = 0; i < n; i ++)
        mpz_init(tmp[i]);

    fft(tmp, n, c, deg);
    for(uint64 i = 0; i < n; i ++) {
        if(mpz_cmp(tmp[i], v[i])) {
            stat = 0;
            break;
        }
    }

    for(uint64 i = 0; i < n; i ++)
        mpz_clear(tmp[i]);
    free(tmp);

    return 1 - stat; //Returns 0 if c and v matches, and 1 otherwise.
}

void coeff_evaluate(mpz_t val, const mpz_t x, const mpz_t *c, const uint64 deg)
{
    mpz_t po, rop, tmp;
    mpz_inits(po, rop, tmp, NULL);

    mpz_set_ui(po, 1);
    mpz_set_ui(rop, 0);
    for(uint64 i = 0; i < deg; i ++) {
        if(i)
            mod_mult(po, po, x);
        mod_mult(tmp, c[i], po);
        mod_add(rop, rop, tmp);
    }

    mpz_set(val, rop);
    mpz_clears(po, rop, tmp, NULL);
}

void polynomial_extrapolate(mpz_t val, const mpz_t x, const mpz_t *xs, const mpz_t *ys, const uint64 n)
{
    mpz_t mult, tmp, inv;
    mpz_inits(mult, tmp, inv, NULL);

    mpz_set_ui(val, 0);
    for(uint64 i = 0; i < n; i ++) {
        mpz_set_ui(mult, 1);
        for(uint64 j = 0; j < n; j ++) {
            if(i != j) {
                mod_sub(tmp, xs[i], xs[j]);
                mod_inv(inv, tmp);
                mod_sub(tmp, x, xs[j]);
                mod_mult(mult, mult, tmp);
                mod_mult(mult, mult, inv);
            }
        }
        mod_mult(tmp, mult, ys[i]);
        mod_add(val, val, tmp);
    }

    mpz_clears(mult, tmp, inv, NULL);
}

void polynomial_extrapolate_N(mpz_t val, const mpz_t x, const mpz_t *ys, const uint64 n)
{
    mpz_t mult, tmp, inv;
    mpz_inits(mult, tmp, inv, NULL);

    mpz_set_ui(val, 0);
    for(uint64 i = 0; i < n; i ++) {
        mpz_set_ui(mult, 1);
        for(uint64 j = 0; j < n; j ++) {
            if(i > j) {
                mod_set_ui(tmp, i - j);
                mod_inv(inv, tmp);
                mod_sub_ui(tmp, x, j);
                mod_mult(mult, mult, tmp);
                mod_mult(mult, mult, inv);
            }
            if(j > i) {
                mod_set_ui(tmp, j - i);
                mod_inv(inv, tmp);
                mod_ui_sub(tmp, j, x);
                mod_mult(mult, mult, tmp);
                mod_mult(mult, mult, inv);
            }
        }
        mod_mult(tmp, mult, ys[i]);
        mod_add(val, val, tmp);
    }

    mpz_clears(mult, tmp, inv, NULL);
}

void fourier_extrapolate(mpz_t val, const mpz_t x, const mpz_t *v, const uint64 n)
{
    mpz_t *xs = (mpz_t*) malloc(sizeof(mpz_t) * n);
    for(uint64 i = 0; i < n; i ++)
        mpz_init_set(xs[i], ROU[i * 4 * N / n]);
        
    polynomial_extrapolate(val, x, xs, v, n);

    for(uint64 i = 0; i < n; i ++)
        mpz_clear(xs[i]);
    free(xs);
}
