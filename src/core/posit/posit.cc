#include <limits.h>
#include "singa/core/posit.h"
#include "posit_helper.h"

using namespace std;

namespace singa {


posit_t posit_zero;
posit_t posit_one;
posit_t posit_min_value;
posit_t posit_max_value;
posit_t posit_nar;

void posit_init() {
    std::cout << "Posit init!" << std::endl;

    posit_zero.ps = DEFAULT_PS;
    posit_zero.es = DEFAULT_ES;    
    SET_POSIT_ZERO(posit_zero.flags);

    posit_nar.ps = DEFAULT_PS;
    posit_nar.es = DEFAULT_ES;        
    SET_POSIT_NAR(posit_nar.flags);
    
    // posit_one = posit_from_int(1);

    int max_rs = DEFAULT_PS - DEFAULT_ES - 1;
    // maximum positive posit  
    posit_max_value.ps = DEFAULT_PS;
    posit_max_value.es = DEFAULT_ES;
    SET_POSIT_SIGN_POSITIVE(posit_max_value.flags);
    posit_max_value.e = ((1 << DEFAULT_ES) - 1);
    posit_max_value.r = max_rs - 2;
    posit_max_value.exp = posit_max_value.r * (1 << DEFAULT_ES) + posit_max_value.e;
    posit_max_value.fs = 0;
    posit_max_value.f = 0;

    // minimum positive posit
    posit_min_value.ps = DEFAULT_PS;
    posit_min_value.es = DEFAULT_ES;
    SET_POSIT_SIGN_POSITIVE(posit_min_value.flags);
    posit_min_value.e = 0;
    posit_min_value.r = 1 - max_rs;
    posit_min_value.exp = posit_min_value.r * (1 << DEFAULT_ES) + posit_min_value.e;
    posit_min_value.fs = 0;
    posit_min_value.f = 0;
}


double posit_to_double(posit_t x) {
    if (IS_POSIT_ZERO(x.flags)) return 0.0;
    if (IS_POSIT_NAR(x.flags)) return NAN;

    // fexp - value of exponent in fractional form
    double fexp = (x.exp < 0) 
        ? 1.0 / (1 << (-x.exp)) 
        : 1.0 * (1 << x.exp);
    int sign = GET_POSIT_SIGN(x.flags);   
    double val = sign * fexp * (1.0 + ((double)x.f)/(1 << x.fs));
    return val;
}

posit_t double_to_posit(double x, uint8_t ps, uint8_t es) {

    // special case - deal fast
    if (x == 0.0)
        return posit_zero;

    // Note: (x == NAN) always return false https://stackoverflow.com/questions/15268864/how-to-compare-two-nan-values-in-c
    if (isnan(x)) {
        return posit_nar;
    }

    double val = ((x<0) * -x) + ((x>=0) * x);
    uint8_t sign_bit = x < 0.0;

    // exp - total exponent: exp = r * 2^es + e
    int exp = (int)log2(val);
    if (exp <= 0 && exp > log2(val)) 
        exp--;
    // e_max - the maximum value of the exponent field (2^es)
    int e_max = (1<<es);
    // r - value of regime field in posit
    int r = exp / e_max;
    // e - value of exponent field in posit
    int e = exp - r * e_max;
    
    // // Correction for the case where exp is neg and is NOT a mult of e_max
    // r = e < 0 ? r - 1 : r;
    // e = exp - r * e_max;
	if (e < 0) {
		r--;
		e = exp - r * e_max;
	}
	if (e > e_max) {
		r++;
		e = exp - r * e_max;
	}


    // The following might not be necessary
    // if (e >= e_max) {
    //     r++;
    //     e = exp - r * e_max;
    // }

//     if (k >= ps - 2) {
// #ifdef PDEBUG
//         printf("Supra-maximal %f!\n", val);
// #endif        
//         if (sign_bit == 1) {
//             return posit_neg(posit_maximum);
//         }
//         else
//             return posit_maximum;
//     }
//     if (k <= -ps + 2) {
// #ifdef PDEBUG
//         printf("Sub-minimal!\n");
// #endif        
//         if (sign == -1)
//             return posit_neg(posit_minimum);
//         else
//             return posit_minimum;
//     }

    // fexp - value of exponent in fractional form
    double fexp = (exp < 0)
        ? 1.0 / (1 << (-exp))
        : 1 << exp;
    // fs - size of fraction
    int fs = (r < 0) 
        ? ps + r - 2 - es   // rs = -(r+1), sign_bit = 1
        : ps - r - 3 - es;  // rs = (r+2),  sign_bit = 1
    int f = 0;
    double f_value = (val / fexp - 1.0);
    if (fs > 0)
        f = f_value * (1 << fs); // f_value = f / (1 << fs)
    else
        fs = 0;

    posit_t ret;
    ret.flags = 0;
    if (sign_bit)
        SET_POSIT_SIGN_NEGATIVE(ret.flags);
    else
        SET_POSIT_SIGN_POSITIVE(ret.flags);
    ret.ps = ps;
    ret.es = es;
    ret.fs = fs;
    ret.r = r;
    ret.e = e;
    ret.exp = exp;
    ret.f = f;

    return ret;
}

posit_t double_to_posit(double x) {
    return double_to_posit(x, DEFAULT_PS, DEFAULT_ES);
}

bool posit_is_equal(posit_t x, posit_t y) {
    return x.ps == y.ps &&
        x.es == y.es &&
        x.fs == y.fs && // optional
        x.flags == y.flags &&
        x.r == y.r &&  // optional
        x.e == y.e &&  // optional
        x.f == y.f &&
        x.exp == y.exp;
}

bool posit_is_less_than(posit_t a, posit_t b) {
    int sign_a = GET_POSIT_SIGN(a.flags);
    int sign_b = GET_POSIT_SIGN(b.flags);
    int both_pos = sign_a == 1 && sign_b == 1;
    int both_neg = sign_a == -1 && sign_b == -1;
    return sign_a < sign_b ||
        (both_neg && a.exp > b.exp) ||
        (both_pos && a.exp < b.exp) ||
        (both_neg && a.exp == b.exp && a.f > b.f) ||
        (both_pos && a.exp == b.exp && a.f < b.f);
}

bool posit_is_less_than_equal(posit_t a, posit_t b) {
    return posit_is_equal(a, b) || posit_is_less_than(a, b);
}

bool posit_is_greater_than(posit_t a, posit_t b) {
    return !posit_is_less_than_equal(a, b);
}

bool posit_is_greater_than_equal(posit_t a, posit_t b) {
    return !posit_is_less_than(a, b);
}

posit_t posit_abs(posit_t x) {
    SET_POSIT_SIGN_POSITIVE(x.flags);
    return x;
}

posit_t posit_neg(posit_t x) {
    TOGGLE_POSIT_SIGN(x.flags);
    return x;
}

posit_t posit_max(posit_t a, posit_t b) {
    return posit_is_greater_than(a,b)
        ? a
        : b;
}

posit_t posit_min(posit_t a, posit_t b) {
    return posit_is_less_than(a,b)
        ? a
        : b;
}

posit_t posit_add(posit_t a, posit_t b) {
    if (IS_POSIT_ZERO(a.flags)) {
        return b;
    } else if (IS_POSIT_ZERO(b.flags)) {
        return a;
    }

    // Full numerator
    long numerator_a = ((long)a.f + ((long)1 << a.fs));
    long numerator_b = ((long)b.f + ((long)1 << b.fs));

    // Diff in fraction size
    int fs_diff = b.fs - a.fs;
    // Base fs is the larger fs
    int fs_base = ((fs_diff >= 0) * b.fs) + ((fs_diff < 0) * a.fs);

    numerator_a = numerator_a << ((fs_diff >= 0) * fs_diff);
    numerator_b = numerator_b << ((fs_diff < 0) * -fs_diff);
    
    // Diff in exponent
    int exp_diff = a.exp - b.exp;
    // Base exp is the smaller exp
    int exp_base = ((exp_diff >= 0) * b.exp) + ((exp_diff < 0) * a.exp);

    numerator_a = numerator_a << ((exp_diff >= 0) * exp_diff);
    numerator_b = numerator_b << ((exp_diff < 0) * -exp_diff);
    
#ifdef PDEBUG
    printf("exp_diff: %d\n", exp_diff);
    printf("fs_diff:  %d\n", fs_diff);
#endif

    // Adjust fraction sign
    numerator_a = numerator_a * GET_POSIT_SIGN(a.flags);
    numerator_b = numerator_b * GET_POSIT_SIGN(b.flags);

    // Result of fraction addition
    long numerator_sum = numerator_a + numerator_b; 
#ifdef PDEBUG
    printf("numerator_a: %lx\n", numerator_a);
    printf("numerator_b: %lx\n", numerator_b);
    printf("numerator_s: %lx\n", numerator_sum);
#endif

    // Prepare result
    posit_t res;
    res.flags = 0;
    int is_neg = numerator_sum < 0;
    res.flags |= (is_neg * 1);
    numerator_sum = numerator_sum * (is_neg * -1 + (1-is_neg));
#ifdef PDEBUG
    printf("numerator_s: %lx\n", numerator_sum);
#endif
    int is_zero = (IS_POSIT_ZERO(a.flags) & IS_POSIT_ZERO(b.flags)) | (numerator_sum == 0);    
    res.flags |= (is_zero * 2);

    // Update exp
    int numerator_exp = (int) log2(numerator_sum);
    long mantissa = numerator_sum - ((long)1 << numerator_exp);
    res.exp = numerator_exp + exp_base - fs_base;
#ifdef PDEBUG
    printf("mantissa:  %lx\n", mantissa);
#endif

	res.e = res.exp % (1 << a.es);
	res.r = res.exp / (1 << a.es);
    
    // Adjust e and r for edge cases when exp is negative
    int adjust_r = ((res.exp < 0) * (res.e != 0) * -1);
    res.r = res.r + adjust_r; // Add 1 to r if exp < 0 and e != 0
    int adjust_e = ((res.exp < 0) * (res.e != 0) * 8); 
    res.e = res.e + adjust_e; // Add 8 to e if exp < 0 and e != 0
    
    // rs = -r + 1 if r < 0 else r + 2
    int rs = (res.r < 0) * (-res.r + 1) + (res.r >= 0) * (res.r + 2);
    // Check whether exp exceeds range
    int max_rs = a.ps - a.es - 1;
    int is_nar = rs > max_rs;
    res.flags |= (is_nar * 4);

    res.fs = a.ps - a.es - rs - 1;

    // Adjust fraction for overflowed exponent
    int adjust_f = res.fs - numerator_exp;
    mantissa = mantissa << ((adjust_f > 0) * adjust_f);
    mantissa = mantissa >> ((adjust_f < 0) * -adjust_f);
    res.f = mantissa;
    res.ps = a.ps;
    res.es = a.es;

#ifdef PDEBUG
    printf("adjust_f:  %d\n", adjust_f);
    printf("res.exp:  %d\n", (int)res.exp);
    printf("res.f:    %d\n", (int)res.f);
    printf("res.fs:   %li\n", ((long)1 << res.fs));    
    printf("res.frac: %lf\n", ((double)res.f / (double)((long)1 << res.fs)));    
#endif

    return res;
}

posit_t posit_sub(posit_t a, posit_t b) {
    return posit_add(a, posit_neg(b));
}

posit_t posit_mul(posit_t a, posit_t b) {
    if (IS_POSIT_ZERO(a.flags) || IS_POSIT_ZERO(b.flags)) {
        return posit_zero;
    } 

    long numerator_a = ((long)a.f + ((long)1 << a.fs));
    long numerator_b = ((long)b.f + ((long)1 << b.fs));
    unsigned long numerator_prod = numerator_a * numerator_b;
    int numerator_exp = log2(numerator_prod);
    long mantissa = numerator_prod - ((long)1 << numerator_exp);

    posit_t res;
    int is_neg = (GET_POSIT_SIGN(a.flags) * GET_POSIT_SIGN(b.flags)) < 0;
    res.flags |= (is_neg * 1);

    res.exp = a.exp + b.exp + numerator_exp - a.fs - b.fs;
	res.e = res.exp % (1 << a.es);
	res.r = res.exp / (1 << a.es);
    
    // Adjust e and r for edge cases when exp is negative
    int adjust_r = ((res.exp < 0) * (res.e != 0) * -1);
    res.r = res.r + adjust_r; // Add 1 to r if exp < 0 and e != 0
    int adjust_e = ((res.exp < 0) * (res.e != 0) * 8); 
    res.e = res.e + adjust_e; // Add 8 to e if exp < 0 and e != 0

    // rs = -r + 1 if r < 0 else r + 2
    int rs = (res.r < 0) * (-res.r + 1) + (res.r >= 0) * (res.r + 2);
    // Check whether exp exceeds range
    int max_rs = a.ps - a.es - 1;
    int is_nar = rs > max_rs;
    res.flags |= (is_nar * 4);

    res.fs = a.ps - a.es - rs - 1;

    // Adjust fraction for overflowed exponent
    int adjust_f = res.fs - numerator_exp;
    mantissa = mantissa << ((adjust_f > 0) * adjust_f);
    mantissa = mantissa >> ((adjust_f < 0) * -adjust_f);
    res.f = mantissa;
    res.ps = a.ps;
    res.es = a.es;

    return res;
}

posit_t posit_div(posit_t a, posit_t b) {
    if (IS_POSIT_ZERO(b.flags)) {
        return posit_nar;
    } 
    if (IS_POSIT_ZERO(a.flags)) {
        return posit_zero;
    }

    long numerator_a = ((long)a.f + ((long)1 << a.fs));
    long numerator_b = ((long)b.f + ((long)1 << b.fs));
    unsigned long numerator_prod = numerator_a * ((long)1 << b.fs) / numerator_b; // this can cause issues
    int numerator_exp = log2(numerator_prod);
    long mantissa = numerator_prod - ((long)1 << numerator_exp);

    posit_t res;
    int is_neg = (GET_POSIT_SIGN(a.flags) * GET_POSIT_SIGN(b.flags)) < 0;
    res.flags |= (is_neg * 1);

    res.exp = a.exp - b.exp + numerator_exp - a.fs;
	res.e = res.exp % (1 << a.es);
	res.r = res.exp / (1 << a.es);
    
    // Adjust e and r for edge cases when exp is negative
    int adjust_r = ((res.exp < 0) * (res.e != 0) * -1);
    res.r = res.r + adjust_r; // Add 1 to r if exp < 0 and e != 0
    int adjust_e = ((res.exp < 0) * (res.e != 0) * 8); 
    res.e = res.e + adjust_e; // Add 8 to e if exp < 0 and e != 0

    // rs = -r + 1 if r < 0 else r + 2
    int rs = (res.r < 0) * (-res.r + 1) + (res.r >= 0) * (res.r + 2);
    // Check whether exp exceeds range
    int max_rs = a.ps - a.es - 1;
    int is_nar = rs > max_rs;
    res.flags |= (is_nar * 4);

    res.fs = a.ps - a.es - rs - 1;

    // Adjust fraction for overflowed exponent
    int adjust_f = res.fs - numerator_exp;
    mantissa = mantissa << ((adjust_f > 0) * adjust_f);
    mantissa = mantissa >> ((adjust_f < 0) * -adjust_f);
    res.f = mantissa;
    res.ps = a.ps;
    res.es = a.es;

    return res;
}


};