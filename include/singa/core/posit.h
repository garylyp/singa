#ifndef SINGA_CORE_POSIT_H_
#define SINGA_CORE_POSIT_H_

#include <iostream>
#include <math.h>
#include <cmath>

namespace singa {

#define DEFAULT_PS  32  // 32-bit posit
#define DEFAULT_ES  3   // 3-bit exponent

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
posit_t posit_sub(posit_t a, posit_t b);
posit_t posit_mul(posit_t a, posit_t b);
posit_t posit_div(posit_t x, posit_t y);




extern posit_t posit_zero;
extern posit_t posit_nar;
extern posit_t posit_max_value;
extern posit_t posit_min_value;

}


#endif  // SINGA_CORE_POSIT_H_