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

    // Adjust fraction of posit with smaller exponent
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
    long mantissa = numerator_sum - ((long)1 << numerator_exp);
    res.exp = numerator_exp + exp_base - fs_base;

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
    printf("res.exp:  %d\n", (int)res.exp);
    printf("res.f:    %d\n", (int)res.f);
    printf("res.fs:   %d\n", ((long)1 << res.fs));    
    printf("res.frac: %lf\n", ((double)res.f / (double)((long)1 << res.fs)));    
#endif

    return res;
}


uint32_t reverse_bits(uint32_t v, int n) {
	uint32_t r = 0;
	int i = 0;

	for (; i < n; i++) {
		r |= ((~(v & 1)) & 0x1) << i;
		v = v >> 1;
	}

	return r;
}

uint32_t nar(int n) {
	return (1 << (n - 1));
}


/**
 * Pack a binary representation of a posit based on relevant fields
 * Inputs:
 * n - posit size
 * es - exponent size
 * special - special posit (0 or NaR)
 * sign - 0 positive, 1 negative
 * k - regime value
 * e - exponent field value
 * fs - fraction size
 * f - fraction value
 * Outputs:
 * posit - the binary representation of the posit
 * err - error code (0 for no error)
 */
void pack_posit(posit_t p, uint32_t* posit, int* err) {
    int n = p.ps;
    int es = p.es;
    int sign = IS_POSIT_NEG(p.flags) ? 1 : 0;
    int special = IS_POSIT_ZERO(p.flags) ? POSIT_ZERO : IS_POSIT_NAR(p.flags) ? POSIT_NAR : POSIT_NOT_SPECIAL;
    int k = p.r;
    int e = p.e;
    int fs = p.fs;
    int f = p.f;

#ifdef PDEBUG
	printf("Pack: n %d, es %d, k %d, e %d, fs %d, f %u\n", n, es, k, e, fs, f);
#endif

	*err = no_err;

	if (special == POSIT_ZERO) {
		*posit = 0;
		return;
	}
	if (special == POSIT_NAR) {
		*posit = nar(n);
		return;
	}

	// x - posit binary value
	uint32_t x = 0;

	int exp_val = (1<<es);
	int bit_index = 0;
	int rs = 0;
	if (k < 0) {
		rs = -k + 1;
		bit_index = n + k - 2;
		x |= (0x1 << bit_index);
	}
	else {
		rs = k + 2;
		bit_index = n - k - 3;
		x |= (((1 << (k + 2)) - 2) << bit_index);
	}

	if (bit_index - es < 0) {
#ifdef PDEBUG
		printf("Error: no space left for the exponent.\n");
#endif
		if (sign)
			x = reverse_bits(x-1, n);
		*posit = x;
		*err = exponent_space;
		return;
	}

	// put the exponent
	x |= (e & (exp_val - 1)) << (bit_index - es);

	// compute final fraction
	int fs2 = n - 1 - rs - es;
	int f2 = f;
	if (fs2 > fs) {
#ifdef PDEBUG
		printf("Difference between given and computed fs: %d %d\n", fs, fs2);
#endif
		f2 = f * (1 << (fs2 - fs));
	}
	else if (fs2 < fs) {
#ifdef PDEBUG
		printf("Difference between given and computed fs: %d %d\n", fs, fs2);
#endif
		f2 = f * (1 << (fs - fs2));
	}

	x |= (f2 & ((1 << fs2) - 1));

	// if the value is negative, take 2's complement
	if (sign)
		x = reverse_bits(x-1, n);

	*posit = x;
}

/**
 * Convert double to posit.
 * Inputs:
 * n - number of posit bits
 * es - number of exponent bits
 * val - double value
 * Outputs:
 * err - error code (0 for no error)
 */
uint32_t encode_posit(int n, int es, double val, int* err) {
	*err = no_err;

	if (n > 8 * sizeof(uint32_t)) {
		*err = posit_space;
#ifdef PDEBUG
		printf("Error: The number of bits n (%d) is larger than the posit storage (%ld).\n", n, 8*sizeof(uint32_t));
#endif
		return 0;
	}

	// x - temporary posit, xx - final posit
	uint32_t x = 0, xx;

	int sign = 0;

	// ival - initial double value (save it)
	double ival = val;

	if (val < 0.0) {
		sign = -1;
		val = -val;
	}
	// special case - deal fast
	if (val == 0.0)
		return 0;
	if (val == NAN)
		return nar(n);	// posit NaR

	// exp - total exponent
	int exp = (int)log2(val);
	if (exp <= 0 && exp > log2(val))
		exp--;
	// exp_val - the maximum value of the exponent field (2^es)
	int exp_val = (1<<es);
	// k - value of regime field in posit
	int k = exp / exp_val;
	// e - value of exponent field in posit
	int e = exp - k * exp_val;
	if (e < 0) {
		k--;
		e = exp - k * exp_val;
	}
	if (e > exp_val) {
		k++;
		e = exp - k * exp_val;
	}
	if (k >= n - 2) {
#ifdef PDEBUG
		printf("Supra-maximal %f!\n", val);
#endif
		*err = supra_max;
		xx = (1 << (n-1)) - 1;
		if (sign == -1) {
			return reverse_bits(xx-1, n);
		}
		else
			return xx;
	}
	if (k <= -n + 2) {
#ifdef PDEBUG
		printf("Sub-minimal!\n");
#endif
		*err = sub_min;
		xx = 1;
		if (sign == -1)
			return reverse_bits(xx-1, n);
		else
			return xx;
	}

	double dexp = (1 << exp);
	if (exp < 0)
		dexp = 1.0 / (1 << (-exp));
	int fs = (k < 0) ? n + k - 2 - es : n - k - 3 - es;
	int frac = 0;
	if (fs > 0)
		frac = (val / dexp - 1.0) * (1 << fs);
	else
		fs = 0;
    posit_t p;
    p.ps = n;
    p.es = es;
    p.flags = 0;
    if (sign) SET_POSIT_SIGN_NEGATIVE(p.flags);
    p.r = k;
    p.e = e;
    p.fs = fs;
    p.f = frac;
	pack_posit(p, &xx, err);

#ifdef PDEBUG
	printf("%x\n", xx);
#endif

	return xx;
}


};