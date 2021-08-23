#ifndef _KERNELS_H_
#define _KERNELS_H_

#include "singa/core/common.h"

namespace singa {

float my_fabs(float x);

float my_sqrt(float x);

float my_exp(float x);

/**
 * Returns the index of the element with the largest absolute value in a vector (single-precision).
 * https://developer.apple.com/documentation/accelerate/1513053-cblas_isamax?language=objc
 */
int my_cblas_isamax(const int __N, const float *__X, const int __incX);

/**
 * Computes the sum of the absolute values of elements in a vector (single-precision).
 * https://developer.apple.com/documentation/accelerate/1513069-cblas_sasum?language=objc
 */
float my_cblas_sasum(const int __N, const float *__X, const int __incX);

/**
 * Computes a constant times a vector plus a vector (single-precision).
 * On return, the contents of vector Y are replaced with the result. The value computed is (alpha * X[i]) + Y[i].
 * https://developer.apple.com/documentation/accelerate/1513188-cblas_saxpy?language=objc
 */
void my_cblas_saxpy(const int __N, const float __alpha, const float *__X, const int __incX, float *__Y, const int __incY);

/**
 * Computes the dot product of two vectors (single-precision).
 * https://developer.apple.com/documentation/accelerate/1513280-cblas_sdot?language=objc
 */
float my_cblas_sdot(const int __N, const float *__X, const int __incX, const float *__Y, const int __incY);

/**
 * Multiplies each element of a vector by a constant (single-precision).
 * https://developer.apple.com/documentation/accelerate/1513354-cblas_sscal?language=objc
 */
void my_cblas_sscal(const int __N, const float __alpha, float *__X, const int __incX);

/**
 * Computes the L2 norm (Euclidian length) of a vector (single precision).
 * https://developer.apple.com/documentation/accelerate/1513250-cblas_snrm2?language=objc
 */
float my_cblas_snrm2(const int __N, const float *__X, const int __incX);

enum MY_CBLAS_ORDER {MyCblasRowMajor=101, MyCblasColMajor=102 };

enum MY_CBLAS_TRANSPOSE {MyCblasNoTrans=111, MyCblasTrans=112, MyCblasConjTrans=113, MyAtlasConj=114};

/**
 * Multiplies a matrix by a vector (single precision).
 * https://developer.apple.com/documentation/accelerate/1513065-cblas_sgemv?language=objc
 */
void my_cblas_sgemv(const enum MY_CBLAS_ORDER __Order, const enum MY_CBLAS_TRANSPOSE __TransA, const int __M, const int __N, const float __alpha, const float *__A, const int __lda, const float *__X, const int __incX, const float __beta, float *__Y, const int __incY);

/**
 * Multiplies two matrices (single-precision).
 * https://developer.apple.com/documentation/accelerate/1513264-cblas_sgemm?language=objc
 */
void my_cblas_sgemm(const enum MY_CBLAS_ORDER __Order, const enum MY_CBLAS_TRANSPOSE __TransA, const enum MY_CBLAS_TRANSPOSE __TransB, const int __M, const int __N, const int __K, const float __alpha, const float *__A, const int __lda, const float *__B, const int __ldb, const float __beta, float *__C, const int __ldc);

} // namespace singa

#endif // _KERNELS_H_
