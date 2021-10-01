/*
 * This file is part of the posit library.
 *
 * Copyright (c) 2019-2020 Dumitrel Loghin.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "conversion.h"

#define EPS 0.0000001

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

char* get_error(int err) {
	switch (err) {
	case wrong_params:
		return "Wrong function parameters";
	case wrong_binary:
		return "Wrong posit binary representation";
	case exponent_space:
		return "No space left for exponent";
	case posit_space:
		return "Not enough space for the posit";
	case supra_max:
		return "Supra-maximal (value greater than maximum psoit)";
	case sub_min:
		return "Sub-minimal (value smaller than minimum posit)";
	case special_posit:
		return "Special posit (0 or NaR) - not an error";
	}
	return "No error";
};

/**
 * Unpack a binary representation of a posit to get relevant fields
 * Inputs:
 * n - posit size
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
void unpack_posit(int n, int es, uint32_t posit, int* special, int* sign, int* k, int* e, int* exp, int* fs, int* f, int* err) {
	*err = no_err;
	*special = POSIT_NOT_SPECIAL;

	// sign
	*sign = 0;
	if (posit & (1<<(n-1))) {
		// negative -> take 2's complement
		*sign = 1;
		posit = reverse_bits(posit, n);
		posit++;
	}

	if (posit == 0) {
		*special = POSIT_ZERO;
		*err = special_posit;
		return;
	}
	if (posit == (1 << (n-1))) {
		*special = POSIT_NAR;
		*err = special_posit;
		return;
	}

	// special case
	if (posit == 1) {
		if (*sign == 0) {
			*k = -n + 2;
		}
		else {
			*k = n - 3;
		}
		*e = 0;
		*fs = 0;
		*f = 0;
		*exp = *k * (1<<es) + *e;
#ifdef PDEBUG
		printf("Unpack(1): %x, n %d, es %d, k %d, e %d, fs %d, f %d\n", posit, n, es, *k, *e, *fs, *f);
#endif
		return;
	}

	// regime
	int regime_value = (posit >> (n-2)) & 0x1;
	// start with the regime field (0-based index => -1, -1 for sign, -1 for first regime field)
	int bit_index = n-3;
	while (bit_index >= 0 && ((posit >> bit_index) & 0x1) == regime_value) bit_index--;

	if (bit_index <= 0) {
#ifdef PDEBUG
		printf("Error: Wrong representation, regime occupies all bits (%08x).\n", posit);
#endif
		*k = (regime_value == 1) ? (n - 3) : (- n + 2);
		*e = 0;
		*exp = *k * (1<<es) + *e;
		*fs = 0;
		*f = 0;
		*err = wrong_binary;
		return;
	}

	// k - value of regime
	*k = (regime_value == 1) ? (n - 3 - bit_index) : (bit_index - n + 2);

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
	printf("Unpack: n %d, es %d, k %d, e %d, fs %d, f %d\n", n, es, *k, *e, *fs, *f);
#endif
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
void pack_posit(int n, int es, int special, int sign, int k, int e, int fs, int f, uint32_t* posit, int* err) {
#ifdef PDEBUG
	printf("Pack: n %d, es %d, k %d, e %d, fs %d, f %d\n", n, es, k, e, fs, f);
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
 * Convert posit to double.
 * Inputs:
 * n - number of bits in the representation
 * es - number of exponent bits
 * posit - posit bits
 */
double decode_posit(int n, int es, uint32_t posit, int* err) {
	int special, sign, k, e, exp, fs, f;
	unpack_posit(n, es, posit, &special, &sign, &k, &e, &exp, &fs, &f, err);

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
	pack_posit(n, es, POSIT_NOT_SPECIAL, sign, k, e, fs, frac, &xx, err);

#ifdef PDEBUG
	printf("%x\n", xx);
#endif

	return xx;
}


/**
 * Transform posit(n1,es1) to posit(n2, es2) where n2 > n1, es2 > es1
 * Inputs:
 * n1 - posit size
 * es1 - exponent size
 * val1 - value of small posit
 * n2 - posit size
 * es2 - exponent size
 * Outputs:
 * err - error code (0 for no error)
 */
uint32_t small_to_big_posit(int n1, int es1, uint32_t val1, int n2, int es2, int* err) {
	*err = no_err;

	if (n1 > n2 || es1 > es2) {
		*err = wrong_params;
		return 0;
	}

	int special1, sign1, k1, e1, exp1, fs1, f1;
	unpack_posit(n1, es1, val1, &special1, &sign1, &k1, &e1, &exp1, &fs1, &f1, err);

	if (special1 == POSIT_ZERO)
		return 0;
	if (special1 == POSIT_NAR)
		return nar(n2);

	if (*err != 0)
		return 0;

	// convert
	uint32_t x = 0;
	int k2 = exp1 / (1 << es2);
	int e2 = exp1 % (1 << es2);
	if (e2 < 0) {
		k2--;
		e2 = exp1 - k2 * (1 << es2);
	}
	int rs2 = 0;
	int fs2 = 0;
	if (k2 < 0) {
		rs2 = -k2 + 1;
	}
	else {
		rs2 = k2 + 2;
	}
	fs2 = n2 - 1 - rs2 - es2;
	int f2 = f1 * (1 << (fs2 - fs1));

	pack_posit(n2, es2, POSIT_NOT_SPECIAL, sign1, k2, e2, fs2, f2, &x, err);

	return x;
}


/**
 * Transform posit(n1,es1) to posit(n2, es2) where n2 < n1, es1 < es2
 * Inputs:
 * n1 - posit size
 * es1 - exponent size
 * val1 - value of small posit
 * n2 - posit size
 * es2 - exponent size
 * Outputs:
 * err - error code (0 for no error)
 */
uint32_t big_to_small_posit(int n1, int es1, uint32_t val1, int n2, int es2, int* err) {
	*err = no_err;

	if (n2 > n1 || es2 > es1) {
		*err = 1;
		return 0;
	}

	int special1, sign1, k1, e1, exp1, fs1, f1;
	unpack_posit(n1, es1, val1, &special1, &sign1, &k1, &e1, &exp1, &fs1, &f1, err);

	if (special1 == POSIT_ZERO)
		return 0;
	if (special1 == POSIT_NAR)
		return nar(n2);

	if (*err != 0)
		return 0;

	// convert
	uint32_t x = 0;
	int k2 = exp1 / (1 << es2);
	int e2 = exp1 % (1 << es2);
	if (e2 < 0) {
		k2--;
		e2 = exp1 - k2 * (1 << es2);
	}
	int rs2 = 0;
	int fs2 = 0;
	if (k2 <= -n2 + 2) {
		if (sign1)
			return reverse_bits(1, n2);
		return 1;
	}
	if (k2 >= n2 - 2) {
		x = (1 << (n2-1)) - 1;
		if (sign1)
			return reverse_bits(x, n2);
		return x;
	}
	if (k2 < 0) {
		rs2 = -k2 + 1;
	}
	else {
		rs2 = k2 + 2;
	}

	//if (rs2 > n2-1-es2) {
	//	*err = 1;
	//	return 0;
	//}

	fs2 = n2 - 1 - rs2 - es2;
	int f2 = f1 / (1 << (fs1 - fs2));

	pack_posit(n2, es2, POSIT_NOT_SPECIAL, sign1, k2, e2, fs2, f2, &x, err);

	return x;
}

