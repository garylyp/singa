#include "posit_helper.h"


namespace singa {

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
 * Unpack a binary representation of a posit to get relevant fields
 * Inputs:
 * ps - posit size
 * es - exponent size
 * posit - the binary representation of the posit
 * Outputs:
 * special - special posit (0 or NaR)
 * sign - 0 positive, 1 negative
 * k - regime value
 * e - exponent field value
 * exp - final exponent value (k*2^es+e)
 * fs - fraction size
 * f - fraction value
 * err - error code (0 for no error)
 */
void unpack_posit(int ps, int es, uint32_t posit, int* special, int* sign, int* k, int* e, int* exp, int* fs, int* f, int* err) {
	*err = no_err;
	*special = POSIT_NOT_SPECIAL;

	// sign
	*sign = 0;
	if (posit & (1<<(ps-1))) {
		// negative -> take 2's complement
		*sign = 1;
		posit = reverse_bits(posit, ps);
		posit++;
	}

	if (posit == 0) {
		*special = POSIT_ZERO;
		*err = special_posit;
		return;
	}
	if (posit == (1 << (ps-1))) {
		*special = POSIT_NAR;
		*err = special_posit;
		return;
	}

	// special case
	if (posit == 1) {
		if (*sign == 0) {
			*k = -ps + 2;
		}
		else {
			*k = ps - 3;
		}
		*e = 0;
		*fs = 0;
		*f = 0;
		*exp = *k * (1<<es) + *e;
#ifdef PDEBUG
		printf("Unpack(1): %x, ps %d, es %d, k %d, e %d, fs %d, f %d\n", posit, ps, es, *k, *e, *fs, *f);
#endif
		return;
	}

	// regime
	int regime_value = (posit >> (ps-2)) & 0x1;
	// start with the regime field (0-based index => -1, -1 for sign, -1 for first regime field)
	int bit_index = ps-3;
	while (bit_index >= 0 && ((posit >> bit_index) & 0x1) == regime_value) bit_index--;

	if (bit_index <= 0) {
#ifdef PDEBUG
		printf("Error: Wrong representation, regime occupies all bits (%08x).\n", posit);
#endif
		*k = (regime_value == 1) ? (ps - 3) : (- ps + 2);
		*e = 0;
		*exp = *k * (1<<es) + *e;
		*fs = 0;
		*f = 0;
		*err = wrong_binary;
		return;
	}

	// k - value of regime
	*k = (regime_value == 1) ? (ps - 3 - bit_index) : (bit_index - ps + 2);

	// e - value of exponent
	*e = ((1 << es) - 1) & (posit >> (bit_index - es));

	// exp - final value of the exponent
	*exp = *k * (1<<es) + *e;

	// fraction size
	bit_index = bit_index - es;
	*fs = bit_index;

	if (*fs < 0) {
#ifdef PDEBUG
		printf("Error: Wrong representation, regime occupies all bits (%08x).\n", posit);
#endif
		*err = wrong_binary;
		return;
	}

	// f - fraction
	*f = ((1 << bit_index) - 1) & posit;

#ifdef PDEBUG
	printf("Unpack: ps %d, es %d, k %d, e %d, fs %d, f %d\n", ps, es, *k, *e, *fs, *f);
#endif
}




/**
 * Pack a binary representation of a posit based on relevant fields
 * Inputs:
 * ps - posit size
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
    int ps = p.ps;
    int es = p.es;
    int sign = IS_POSIT_NEG(p.flags) ? 1 : 0;
    int special = IS_POSIT_ZERO(p.flags) ? POSIT_ZERO : IS_POSIT_NAR(p.flags) ? POSIT_NAR : POSIT_NOT_SPECIAL;
    int k = p.r;
    int e = p.e;
    int fs = p.fs;
    int f = p.f;

#ifdef PDEBUG
	printf("Pack: ps %d, es %d, k %d, e %d, fs %d, f %u\n", ps, es, k, e, fs, f);
#endif

	*err = no_err;

	if (special == POSIT_ZERO) {
		*posit = 0;
		return;
	}
	if (special == POSIT_NAR) {
		*posit = nar(ps);
		return;
	}

	// x - posit binary value
	uint32_t x = 0;

	int exp_val = (1<<es);
	int bit_index = 0;
	int rs = 0;
	if (k < 0) {
		rs = -k + 1;
		bit_index = ps + k - 2;
		x |= (0x1 << bit_index);
	}
	else {
		rs = k + 2;
		bit_index = ps - k - 3;
		x |= (((1 << (k + 2)) - 2) << bit_index);
	}

	if (bit_index - es < 0) {
#ifdef PDEBUG
		printf("Error: no space left for the exponent.\n");
#endif
		if (sign)
			x = reverse_bits(x-1, ps);
		*posit = x;
		*err = exponent_space;
		return;
	}

	// put the exponent
	x |= (e & (exp_val - 1)) << (bit_index - es);

	// compute final fraction
	int fs2 = ps - 1 - rs - es;
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
		x = reverse_bits(x-1, ps);

	*posit = x;
}

/**
 * Convert double to posit.
 * Inputs:
 * ps - number of posit bits
 * es - number of exponent bits
 * val - double value
 * Outputs:
 * err - error code (0 for no error)
 */
uint32_t encode_posit(int ps, int es, double val, int* err) {
	*err = no_err;

	if (ps > 8 * sizeof(uint32_t)) {
		*err = posit_space;
#ifdef PDEBUG
		printf("Error: The number of bits n (%d) is larger than the posit storage (%ld).\n", ps, 8*sizeof(uint32_t));
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
		return nar(ps);	// posit NaR

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
	if (k >= ps - 2) {
#ifdef PDEBUG
		printf("Supra-maximal %f!\n", val);
#endif
		*err = supra_max;
		xx = (1 << (ps-1)) - 1;
		if (sign == -1) {
			return reverse_bits(xx-1, ps);
		}
		else
			return xx;
	}
	if (k <= -ps + 2) {
#ifdef PDEBUG
		printf("Sub-minimal!\n");
#endif
		*err = sub_min;
		xx = 1;
		if (sign == -1)
			return reverse_bits(xx-1, ps);
		else
			return xx;
	}

	double dexp = (1 << exp);
	if (exp < 0)
		dexp = 1.0 / (1 << (-exp));
	int fs = (k < 0) ? ps + k - 2 - es : ps - k - 3 - es;
	int frac = 0;
	if (fs > 0)
		frac = (val / dexp - 1.0) * (1 << fs);
	else
		fs = 0;
    posit_t p;
    p.ps = ps;
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


double decode_posit(int ps, int es, uint32_t posit, int* err) {
	int special, sign, k, e, exp, fs, f;
	unpack_posit(ps, es, posit, &special, &sign, &k, &e, &exp, &fs, &f, err);

	if (special == POSIT_ZERO)
		return 0.0;
	if (special == POSIT_NAR)
		return NAN;

	if (*err != no_err)
		return NAN;

	if (sign == 0)
		sign = 1;
	else
		sign = -1;

	double fexp = (exp < 0) ? 1.0 / (1 << (-1 * exp)) : 1.0 * (1 << exp);

	// x - converted value
	double x = sign * fexp * (1.0 + ((double)f)/(1 << fs));

#ifdef PDEBUG
	printf("k %i, e %i, exp %i, f %x\n", k, e, exp, f);
	printf("%15.14lf\n", x);
#endif

	return x;
}

void compare_posit(uint32_t a, uint32_t b, int ps, int es) {
    int a_special, a_sign, a_r, a_e, a_exp, a_fs, a_f, a_err;
    int b_special, b_sign, b_r, b_e, b_exp, b_fs, b_f, b_err;
	printf("Expected:\n");
    unpack_posit(ps, es, a, &a_special, &a_sign, &a_r, &a_e, &a_exp, &a_fs, &a_f, &a_err);
	printf("Actual:\n");
    unpack_posit(ps, es, b, &b_special, &b_sign, &b_r, &b_e, &b_exp, &b_fs, &b_f, &b_err);
    if (a_err) {
        printf("error unpacking a: %d", a_err);
    }
    if (b_err) {
        printf("error unpacking b: %d", b_err);
    }
    printf("a_sign: %d, b_sign: %d\n", a_sign, b_sign);
    printf("a_r: %d, b_r: %d\n", a_r, b_r);
    printf("a_e: %d, b_e: %d\n", a_e, b_e);
    printf("a_exp: %d, b_exp: %d\n", a_exp, b_exp);
    printf("a_fs: %d, b_fs: %d\n", a_fs, b_fs);
    printf("a_f: %x, b_f: %x\n", a_f, b_f);
    printf("a_special: %d, b_special: %d\n", a_special, b_special);
    printf("special: POSIT_ZERO = 0, POSIT_NAR = 1, POSIT_NOT_SPECIAL = 2\n");
}

void unpack_exp(int exp, int es, int* r, int* e) {
	int e_val, r_val;
	e_val = exp % (1 << es);
	r_val = exp / (1 << es);
    
    // Adjust e and r for edge cases when exp is negative
    int adjust_r = ((exp < 0) * (e_val != 0) * -1);
    r_val = r_val + adjust_r; // Add 1 to r if exp < 0 and e != 0
    int adjust_e = ((exp < 0) * (e_val != 0) * 8); 
    e_val = e_val + adjust_e; // Add 8 to e if exp < 0 and e != 0

    *r = r_val;
	*e = e_val;
}


}