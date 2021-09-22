#include "singa/core/posit.h"
#include <limits.h>

using namespace std;

namespace singa {

#define IS_POSIT_NEG(flags) ((flags & 1) == 1)
#define IS_POSIT_POS(flags) ((flags & 1) == 0)
#define SET_POSIT_SIGN_POSITIVE(flags) flags &= 0xFE
#define SET_POSIT_SIGN_NEGATIVE(flags) flags |= 1
#define TOGGLE_POSIT_SIGN(flags) flags ^= 1
#define GET_POSIT_SIGN(flags) (((flags & 1) == 0) ? 1 : -1)

#define IS_POSIT_ZERO(flags) ((flags & 2) != 0)
#define SET_POSIT_ZERO(flags) flags |= 2
#define SET_POSIT_NONZERO(flags) flags &= 0xFD
#define TOGGLE_POSIT_ZERO(flags) flags ^= 2

#define IS_POSIT_NAR(flags) ((flags & 4) != 0)
#define SET_POSIT_NAR(flags) flags |= 4
#define TOGGLE_POSIT_NAR(flags) flags ^= 4

posit_t posit_zero;
posit_t posit_one;
posit_t posit_minimum;
posit_t posit_maximum;
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

    // // maximum positive posit  
    // posit_maximum.ps = DEFAULT_PS;
    // posit_maximum.es = DEFAULT_ES;
    // SET_POSIT_SIGN_POSITIVE(posit_maximum.flags);
    // posit_maximum.e = ((1 << DEFAULT_ES) - 1);
    // posit_maximum.k = DEFAULT_PS - DEFAULT_ES - 3;
    // posit_maximum.exp = posit_maximum.k * (1 << DEFAULT_ES) + posit_maximum.e;
    // posit_maximum.fs = 0;
    // posit_maximum.f = 0;

    // // minimum positive posit
    // posit_minimum.ps = DEFAULT_PS;
    // posit_minimum.es = DEFAULT_ES;
    // SET_POSIT_SIGN_POSITIVE(posit_minimum.flags);
    // posit_minimum.e = 0;
    // posit_minimum.k = 2 + DEFAULT_ES - DEFAULT_PS;
    // posit_minimum.exp = posit_maximum.k * (1 << DEFAULT_ES) + posit_maximum.e;
    // posit_minimum.fs = 0;
    // posit_minimum.fraction = 0;
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

    double val = abs(x);
    uint8_t sign_bit = x < 0.0;

    // exp - total exponent: exp = r * 2^es + e
    int exp = (int)log2(val);
    if (exp <= 0 && exp > log2(val)) // What is this step?
        exp--;
    // e_max - the maximum value of the exponent field (2^es)
    int e_max = (1<<es);
    // r - value of regime field in posit
    int r = exp / e_max;
    // e - value of exponent field in posit
    int e = exp - r * e_max;
    
    // Correction for the case where exp is neg and is NOT a mult of e_max
    r = e < 0 ? r - 1 : r;
    e = exp - r * e_max;

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
    // // if diff in exp > (ps - es - 1 - 2), just return the one with larger exponent
    // int max_fs_diff = a.ps - 2 - 1 - a.es;
    if (IS_POSIT_ZERO(a.flags)) {
        return b;
    } else if (IS_POSIT_ZERO(b.flags)) {
        return a;
    }
    
    // Diff in exponent
    int exp_diff = a.exp - b.exp;

    // We make 'a' the one with the higher exponent
    posit_t c = a;
    if (exp_diff < 0) {
        a = b;
        b = c;
        exp_diff = -exp_diff;
    } 

    // Adjust to common exponent (b's exp) by raising a's fraction   
    long numerator_a = ((long)a.f + ((long)1 << a.fs)) << exp_diff;
    long numerator_b = (b.f + (1 << b.fs));

    // Adjust to common fraction denominator by lowering a's fraction
    int fs_diff = b.fs - a.fs;
    numerator_a = numerator_a << fs_diff;

    // Adjust fraction sign
    numerator_a = numerator_a * GET_POSIT_SIGN(a.flags);
    numerator_b = numerator_b * GET_POSIT_SIGN(b.flags);

    // Add
    long numerator_sum = numerator_a + numerator_b; // new f, but unadjusted and might be too large / negative

    // Prepare result
    posit_t res;
    res.flags = 0;
    int is_neg = numerator_sum < 0;
    
    res.flags |= (is_neg * 1);
    numerator_sum = numerator_sum * (is_neg * -1 + (1-is_neg));

    int is_zero = (IS_POSIT_ZERO(a.flags) & IS_POSIT_ZERO(b.flags)) | (numerator_sum == 0);    
    res.flags |= (is_zero * 2);

    // Update exp
    int numerator_exp = log2(numerator_sum);
    // over/under-flowed exponent to be shifted from numerator to exp
    int overflowed_exp = numerator_exp - b.fs; 
    res.exp = b.exp + overflowed_exp; 
#ifdef PDEBUG
    cout << "res.exp: " << (int)res.exp << endl; 
#endif
	res.e = res.exp % (1 << a.es);
	res.r = res.exp / (1 << a.es);
    
    // Adjust e and r for edge cases when exp is negative
    int adjust_r = ((res.e == 0) * (res.exp < 0)); 
    res.r = res.r + adjust_r;
    int adjust_e = ((res.e != 0) * (res.exp < 0) * 8); 
    res.e = res.e + adjust_e;
    int rs = (res.r < 0) * (-res.r + 1) + (res.r >= 0) * (res.r + 2);
#ifdef PDEBUG
    cout << "res.r: " << (int)res.r << endl;
    cout << "res.e: " << (int)res.e << endl;
#endif
    // Check whether exp exceeds range
    int max_rs = a.ps - a.es - 1;
    int is_nar = rs > max_rs;
    res.flags |= (is_nar * 4);

    // Adjust fraction for overflowed exponent
    numerator_sum = numerator_sum >> ((overflowed_exp > 0) * overflowed_exp);
    numerator_sum = numerator_sum << ((overflowed_exp <= 0) * -overflowed_exp);
    numerator_sum = numerator_sum - (1 << b.fs); // minus 1 from the fraction (minus denominator off numerator)
    res.ps = a.ps;
    res.es = a.es;
    res.fs = res.ps - res.es - rs;

#ifdef PDEBUG
    cout << "numer: " << numerator_sum << endl;
    cout << "denom: " << (1<<b.fs) << endl;
#endif
    
    // Might unncessarily truncate some items, need more tests
    numerator_sum = numerator_sum << ((res.fs > b.fs) * (res.fs - b.fs));
    numerator_sum = numerator_sum >> ((res.fs <= b.fs) * (b.fs - res.fs));

    res.f = numerator_sum;
#ifdef PDEBUG
    cout << "numer: " << numerator_sum << endl;
    cout << "denom: " << (1<<res.fs) << endl;
#endif

    return res;
}



uint32_t int_log2 (uint32_t val) {
    if (val == 0) return UINT_MAX;
    if (val == 1) return 0;
    uint32_t ret = 0;
    while (val > 1) {
        val >>= 1;
        ret++;
    }
    return ret;
}
};