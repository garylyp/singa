#ifndef SINGA_CORE_POSIT_H_
#define SINGA_CORE_POSIT_H_

#include <iostream>
#include <math.h>
#include <cmath>

namespace singa {

#define DEFAULT_PS  32  // 32-bit posit
#define DEFAULT_ES  3   // 3-bit exponent
struct posit_t;

struct posit_t {

  // ...[IS_NAR][IS_ZERO][IS_NEG]
  uint8_t flags;
  uint8_t ps;
  uint8_t es;
  uint8_t fs;
  int8_t r;
  uint8_t e;
  uint32_t f;
  int8_t exp;
  
  posit_t() {};
  // posit_t(int x) {};
  // posit_t(double x) {};
  // posit_t operator/ (posit_t y) const { return posit_div(*this, y); };
  // bool operator != (posit_t y) const { return !posit_is_equal(*this, y); };
  // bool operator == (posit_t y) const { return posit_is_equal(*this, y); };
  // friend std::ostream& operator<< (std::ostream& out, const posit_t& x) { out << "posit"; return out; };
}; 


void posit_init();
// posit_t posit_from_int(int x);
double posit_to_double(posit_t x);
posit_t double_to_posit(double x, uint8_t ps, uint8_t es);
posit_t double_to_posit(double x);



bool posit_is_equal(posit_t x, posit_t y);
bool posit_is_less_than(posit_t a, posit_t b);
bool posit_is_less_than_equal(posit_t a, posit_t b);
bool posit_is_greater_than(posit_t a, posit_t b);
bool posit_is_greater_than_equal(posit_t a, posit_t b);

posit_t posit_abs(posit_t x);
posit_t posit_neg(posit_t x);
// posit_t posit_sign(posit_t x);

posit_t posit_max(posit_t a, posit_t b);
posit_t posit_min(posit_t a, posit_t b);
posit_t posit_add(posit_t a, posit_t b);
// posit_t posit_sub(posit_t a, posit_t b);
// posit_t posit_mul(posit_t a, posit_t b);
// posit_t posit_div(posit_t x, posit_t y);


enum err_code_t {
	no_err = 0,
	wrong_params = 1,
	wrong_binary = 2,
	exponent_space = 3,
	posit_space = 4,
	supra_max = 5,
	sub_min = 6,
	special_posit = 7
};

enum special_posit_t {
	POSIT_ZERO = 0,
	POSIT_NAR = 1,
	POSIT_NOT_SPECIAL = 2
};
uint32_t encode_posit(int n, int es, double val, int* err);
void pack_posit(posit_t p, uint32_t* posit, int* err);

extern posit_t posit_zero;
extern posit_t posit_nar;
}


#endif  // SINGA_CORE_POSIT_H_