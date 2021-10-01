#ifndef SINGA_CORE_POSIT_H_
#define SINGA_CORE_POSIT_H_

#include <iostream>

namespace singa {

enum MY_CBLAS_ORDER {MyCblasRowMajor=101, MyCblasColMajor=102 };

enum MY_CBLAS_TRANSPOSE {MyCblasNoTrans=111, MyCblasTrans=112, MyCblasConjTrans=113, MyAtlasConj=114};

#define DEFAULT_PS  32
#define DEFAULT_ES  3

struct posit_t;

void posit_init();

int posit_to_int(posit_t x);
posit_t posit_from_int(int x);
float posit_to_float(posit_t x);
posit_t posit_from_float(float x);

// some workarounds
float to_float(posit_t x);
float to_float(int x);

bool posit_is_equal(posit_t x, posit_t y);
bool posit_is_smaller(posit_t a, posit_t b);
bool posit_is_smaller_equal(posit_t a, posit_t b);
bool posit_is_bigger(posit_t a, posit_t b);
bool posit_is_bigger_equal(posit_t a, posit_t b);

posit_t posit_abs(posit_t x);
posit_t posit_neg(posit_t x);
posit_t posit_sign(posit_t x);
posit_t posit_max(posit_t a, posit_t b);
posit_t posit_min(posit_t a, posit_t b);
posit_t posit_add(posit_t a, posit_t b);
posit_t posit_sub(posit_t a, posit_t b);
posit_t posit_mul(posit_t a, posit_t b);
posit_t posit_div(posit_t x, posit_t y);
posit_t posit_sqrt(posit_t x);
posit_t posit_exp(posit_t x);
posit_t posit_log(posit_t x);
posit_t posit_pow(posit_t y, posit_t x);
posit_t posit_tanh(posit_t x);

// int posit_cblas_isamax(const int __N, const posit_t *__X, const int __incX);
// posit_t posit_cblas_sasum(const int __N, const posit_t *__X, const int __incX);
// void posit_cblas_saxpy(const int __N, const posit_t __alpha, const posit_t *__X, const int __incX, posit_t *__Y, const int __incY);
// posit_t posit_cblas_sdot(const int __N, const posit_t *__X, const int __incX, const posit_t *__Y, const int __incY);
// void posit_cblas_sscal(const int __N, const posit_t __alpha, posit_t *__X, const int __incX);
// posit_t posit_cblas_snrm2(const int __N, const posit_t *__X, const int __incX);
// void posit_cblas_sgemv(const enum MY_CBLAS_ORDER __Order,
// 		const enum MY_CBLAS_TRANSPOSE __TransA,
// 		const int __M, const int __N,
// 		const posit_t __alpha,
// 		const posit_t *__A, const int __lda,
// 		const posit_t *__X, const int __incX,
// 		const posit_t __beta,
// 		posit_t *__Y, const int __incY);
// void posit_cblas_sgemm(const enum MY_CBLAS_ORDER __Order,
// 		const enum MY_CBLAS_TRANSPOSE __TransA,
// 		const enum MY_CBLAS_TRANSPOSE __TransB,
// 		const int __M, const int __N, const int __K,
// 		const posit_t __alpha,
// 		const posit_t *__A, const int __lda,
// 		const posit_t *__B, const int __ldb,
// 		const posit_t __beta,
// 		posit_t *__C, const int __ldc);


struct posit_t {
  uint8_t flags;
  uint8_t ps;
  uint8_t es;
  uint8_t fs;
  int8_t k;
  uint8_t e;
  int8_t exp;
  uint32_t fraction;
  
  posit_t() {};
  posit_t(int x) {};
  posit_t operator/ (posit_t y) const { return posit_div(*this, y); };
  bool operator != (posit_t y) const { return !posit_is_equal(*this, y); };
  bool operator == (posit_t y) const { return posit_is_equal(*this, y); };
  friend std::ostream& operator<< (std::ostream& out, const posit_t& x) { out << "posit"; return out; };
}; 

// #define posit_t(X) posit_from_int(X)

extern posit_t posit_zero;
extern posit_t posit_one;
extern posit_t posit_minimum;
extern posit_t posit_maximum;

}

#endif  // SINGA_CORE_POSIT_H_