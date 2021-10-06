
#include <math.h>
#include "posit2.h"

namespace singa {

#define IS_POSIT_ZERO(flags) ((flags & 2) != 0)
#define IS_POSIT_NAR(flags) ((flags & 4) != 0)

#define SET_POSIT_ZERO(flags) flags |= 2
#define SET_POSIT_NAR(flags) flags |= 4

#define GET_POSIT_SIGN(flags) (((flags & 1) == 0) ? 1 : -1)
#define SET_POSIT_SIGN_POSITIVE(flags) flags &= 0xFE
#define SET_POSIT_SIGN_NEGATIVE(flags) flags |= 1

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
    
	posit_one = posit_from_int(1);

	// maximum positive posit  
	posit_maximum.ps = DEFAULT_PS;
	posit_maximum.es = DEFAULT_ES;
	SET_POSIT_SIGN_POSITIVE(posit_maximum.flags);
	posit_maximum.e = ((1 << DEFAULT_ES) - 1);
	posit_maximum.k = DEFAULT_PS - DEFAULT_ES - 3;
	posit_maximum.exp = posit_maximum.k * (1 << DEFAULT_ES) + posit_maximum.e;
	posit_maximum.fs = 0;
	posit_maximum.fraction = 0;

	// minimum positive posit
	posit_minimum.ps = DEFAULT_PS;
	posit_minimum.es = DEFAULT_ES;
	SET_POSIT_SIGN_POSITIVE(posit_minimum.flags);
	posit_minimum.e = 0;
	posit_minimum.k = 2 + DEFAULT_ES - DEFAULT_PS;
	posit_minimum.exp = posit_maximum.k * (1 << DEFAULT_ES) + posit_maximum.e;
	posit_minimum.fs = 0;
	posit_minimum.fraction = 0;
}

bool posit_is_smaller(posit_t a, posit_t b) {
    int sign_a = GET_POSIT_SIGN(a.flags);
    int sign_b = GET_POSIT_SIGN(b.flags);
    if (sign_a == -1 && sign_b == 1)
        return true;
    if (sign_a == 1 && sign_b == -1)
        return false;
    if (a.exp < b.exp)
        return true;
    if (a.exp > b.exp)
        return false;
    if (a.fraction < b.fraction)
        return true;    
    return false;
}

bool posit_is_equal(posit_t a, posit_t b) {
    return (a.ps == b.ps) && (a.es == b.es) && (a.flags == b.flags) && (a.exp == b.exp) && (a.fraction == b.fraction);
}

bool posit_is_smaller_equal(posit_t a, posit_t b) {
    return posit_is_equal(a, b) || posit_is_smaller(a, b);
}

bool posit_is_bigger(posit_t a, posit_t b) {
    return posit_is_smaller(b, a);
}

bool posit_is_bigger_equal(posit_t a, posit_t b) {
    return posit_is_equal(a, b) || posit_is_smaller(b, a);
}

double posit_to_double(posit_t x) {
    if (IS_POSIT_ZERO(x.flags))
		return 0.0;
	if (IS_POSIT_NAR(x.flags))
		return NAN;

    int sign = GET_POSIT_SIGN(x.flags);   
	double fexp = (x.exp < 0) ? 1.0 / (1 << (-1 * x.exp)) : 1.0 * (1 << x.exp);
	double val = sign * fexp * (1.0 + ((double)x.fraction)/(1 << x.fs));

    return val;
}

posit_t posit_from_double(double val, uint8_t ps, uint8_t es) {
    int sign = 0;

	// ival - initial double value (save it)
	double ival = val;

	if (val < 0.0) {
		sign = -1;
		val = -val;
	}
	// special case - deal fast
	if (val == 0.0)
		return posit_zero;
	if (val == NAN)
		return posit_nar;

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
		if (sign == -1) {
			return posit_neg(posit_maximum);
		}
		else
			return posit_maximum;
	}
	if (k <= -ps + 2) {
#ifdef PDEBUG
		printf("Sub-minimal!\n");
#endif		
		if (sign == -1)
			return posit_neg(posit_minimum);
		else
			return posit_minimum;
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

    posit_t ret;
    if (sign == -1)
        SET_POSIT_SIGN_NEGATIVE(ret.flags);
    else
        SET_POSIT_SIGN_POSITIVE(ret.flags);
    ret.ps = ps;
    ret.es = es;
    ret.fs = fs;
    ret.k = k;
    ret.e = e;
    ret.exp = exp;
    ret.fraction = frac;

    return ret;
}

float posit_to_float(posit_t x) {
    return (float)posit_to_double(x);
}

float to_float(posit_t x) {
    return posit_to_float(x);
}

float to_float(int x) {
    return static_cast<float>(x);
}

posit_t posit_from_float(float x) {    
    return posit_from_double((double)x, DEFAULT_PS, DEFAULT_ES);
}

posit_t posit_from_float(float x, uint8_t ps, uint8_t es) {    
    return posit_from_double((double)x, ps, es);
}

posit_t posit_from_int(int x) {
   return posit_from_double(1.0 * x, DEFAULT_PS, DEFAULT_ES);
}

posit_t posit_abs(posit_t x) {
    if (GET_POSIT_SIGN(x.flags) == -1)
        return posit_neg(x);
    return x;
}

posit_t posit_neg(posit_t x) {
	if (GET_POSIT_SIGN(x.flags) == 1)
		SET_POSIT_SIGN_NEGATIVE(x.flags);
	else
		SET_POSIT_SIGN_POSITIVE(x.flags);
    return x;
}

posit_t posit_sign(posit_t x) {
    return (GET_POSIT_SIGN(x.flags) == 1) ? posit_one : posit_neg(posit_one);
}

posit_t posit_max(posit_t a, posit_t b) {
    if (posit_is_bigger(a, b))
        return a;
    return b;
}

posit_t posit_min(posit_t a, posit_t b) {
    if (posit_is_smaller(a, b))
        return a;
    return b;
}

// TODO
posit_t posit_add(posit_t a, posit_t b) {
	double da = posit_to_double(a);
	double db = posit_to_double(b);
	return posit_from_double(da + db, DEFAULT_PS, DEFAULT_ES);
}
// TODO
posit_t xxx_posit_add(posit_t a, posit_t b) {
    if (IS_POSIT_NAR(a.flags) || IS_POSIT_NAR(b.flags))
        return posit_nar;
    if (IS_POSIT_ZERO(a.flags))
        return b;
    if (IS_POSIT_ZERO(b.flags)) {
        return a;
	}
        
	int sign_a = GET_POSIT_SIGN(a.flags);
	int sign_b = GET_POSIT_SIGN(b.flags);
	double f1 = 1.0 + (double)a.fraction / (1 << a.fs);
	double f2 = 1.0 + (double)b.fraction / (1 << b.fs);
	double f3;

	posit_t ret;	
	ret.ps = a.ps;
	ret.es = a.es;

	if (a.exp < b.exp) {
			f1 = f1 / (1 << (b.exp - a.exp));
			ret.exp = b.exp;
		}
		else {
			if (a.exp > b.exp) {
				f2 = f2 / (1 << (a.exp - b.exp));				
			}
			ret.exp = a.exp;			
		}

	if (sign_a == sign_b) {
		if (sign_a == 1)
			SET_POSIT_SIGN_POSITIVE(a.flags);
		else
			SET_POSIT_SIGN_NEGATIVE(b.flags);		
		f3 = f1 + f2;
	}
	else {
		if (sign_a == 1)
			f3 = fabs(f1 - f2);
		else
			f3 = fabs(f2 - f1);
	}

	if (f3 >= 2.0) {
		f3 = f3 / 2.0;
		ret.exp += 1;
	}
	if (f3 < 1.0) {
		f3 = f3 * 2.0;
		ret.exp -= 1;
	}

	ret.e = ret.exp % (1 << a.es);
	ret.k = ret.exp / (1 << a.es);
	if (ret.k < 1 + a.es - a.ps)
		return posit_minimum;
	if (ret.k > a.ps - a.es - 1)
		return posit_maximum;
	int rs = 1 + (ret.k < 0) ? -ret.k : ret.k + 1;	

	ret.fs = ret.ps - ret.es - rs;
	ret.fraction = (uint32_t)((f3 - 1.0) * (1 << ret.fs));

	return ret;
}

// TODO
posit_t posit_sub(posit_t a, posit_t b) {
	double da = posit_to_double(a);
	double db = posit_to_double(b);
	return posit_from_double(da - db, DEFAULT_PS, DEFAULT_ES);
}
// TODO
posit_t xxx_posit_sub(posit_t a, posit_t b) {    
    return posit_add(a, posit_neg(b));
}

// TODO
posit_t posit_mul(posit_t a, posit_t b) {
	double da = posit_to_double(a);
	double db = posit_to_double(b);
	return posit_from_double(da * db, DEFAULT_PS, DEFAULT_ES);
}
// TODO
posit_t xxx_posit_mul(posit_t a, posit_t b) {
    if (IS_POSIT_NAR(a.flags) || IS_POSIT_NAR(b.flags))
        return posit_nar;
    if (IS_POSIT_ZERO(a.flags) || IS_POSIT_ZERO(b.flags))
        return posit_zero;
    if (posit_is_equal(a, posit_one))
        return b;
    if (posit_is_equal(b, posit_one))
        return a;
    int sign = GET_POSIT_SIGN(a.flags) * GET_POSIT_SIGN(b.flags);
    int exp = a.exp + b.exp;
    
	double f1 = 1.0 + (double)a.fraction / (1 << a.fs);
	double f2 = 1.0 + (double)b.fraction / (1 << b.fs);
	double f3 = f1 * f2;
	if (f3 >= 2.0) {
		f3 = f3 / 2.0;
		exp += 1;
	}

	posit_t ret;	
	if (sign == -1)
		SET_POSIT_SIGN_NEGATIVE(ret.flags);
	else
		SET_POSIT_SIGN_POSITIVE(ret.flags);
	ret.ps = a.ps;
	ret.es = a.es;
	ret.e = exp % (1 << a.es);
	ret.k = exp / (1 << a.es);
	if (ret.k < 1 + a.es - a.ps)
		return posit_minimum;
	if (ret.k > a.ps - a.es - 1)
		return posit_maximum;
	int rs = 1 + (ret.k < 0) ? -ret.k : ret.k + 1;
	ret.exp = exp;

	ret.fs = ret.ps - ret.es - rs;
	ret.fraction = (uint32_t)((f3 - 1.0) * (1 << ret.fs));

	return ret;
}

// TODO
posit_t posit_div(posit_t a, posit_t b) {
	double da = posit_to_double(a);
	double db = posit_to_double(b);
	return posit_from_double(da / db, DEFAULT_PS, DEFAULT_ES);
}
// TODO
posit_t xxx_posit_div(posit_t a, posit_t b) {
    if (IS_POSIT_ZERO(b.flags))
        return posit_nar;
    if (IS_POSIT_ZERO(a.flags))
        return posit_zero;
    if (posit_is_equal(b, posit_one))
        return a;

    int sign = GET_POSIT_SIGN(a.flags) * GET_POSIT_SIGN(b.flags);
    int exp = a.exp - b.exp;

	double f1 = 1.0 + (double)a.fraction / (1 << a.fs);
	double f2 = 1.0 + (double)b.fraction / (1 << b.fs);
	double f3 = f1 / f2;
	if (f3 < 1.0) {
		f3 = f3 * 2.0;
		exp -= 1;
	}

	posit_t ret;	
	if (sign == -1)
		SET_POSIT_SIGN_NEGATIVE(ret.flags);
	else
		SET_POSIT_SIGN_POSITIVE(ret.flags);
	ret.ps = a.ps;
	ret.es = a.es;
	ret.e = exp % (1 << a.es);
	ret.k = exp / (1 << a.es);
	if (ret.k < 2 + a.es - a.ps)
		return posit_minimum;
	if (ret.k > a.ps - a.es - 2)
		return posit_maximum;
	int rs = 1 + (ret.k < 0) ? -ret.k : ret.k + 1;
	ret.exp = exp;

	ret.fs = ret.ps - ret.es - rs;
	ret.fraction = (uint32_t)((f3 - 1.0) * (1 << ret.fs));

	return ret;
}

posit_t posit_sqrt(posit_t x) {
    double dx = posit_to_double(x);    
    return posit_from_double(sqrt(dx), x.ps, x.es);
}

posit_t posit_exp(posit_t x) {
    double dx = posit_to_double(x);
    return posit_from_double(exp(dx), x.ps, x.es);
}

posit_t posit_log(posit_t x) {
    double dx = posit_to_double(x);
    return posit_from_double(log(dx), x.ps, x.es);
}

posit_t posit_pow(posit_t x, posit_t y) {
    double dx = posit_to_double(x);
    double dy = posit_to_double(y);
    return posit_from_double(pow(dx, dy), x.ps, x.es);
}

posit_t posit_tanh(posit_t x) {
    double dx = posit_to_double(x);
    return posit_from_double(tanh(dx), x.ps, x.es);
}
}