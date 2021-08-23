#include "kernels.h"

// Constants
namespace singa {

float const_float_zero = 0.f;
float const_float_one = 1.f;
float const_float_minus_one = -1.f;
float const_float_two = 2.f;
float const_float_minus_two = -2.f;
float const_float_zero_five = 0.5f;
float const_float_minus_zero_five = -0.5f;
float const_float_one_five = 1.5f;
float const_float_minus_one_five = -1.5f;
float const_float_1e_6f = 1e-6f;

void init_const_float() {
	const_float_zero = 0.f;
	const_float_one = 1.f;
	const_float_minus_one = -1.f;
	const_float_two = 2.f;
	const_float_minus_two = -2.f;
	const_float_zero_five = 0.5f;
	const_float_minus_zero_five = -0.5f;
	const_float_one_five = 1.5f;
	const_float_minus_one_five = -1.5f;
	const_float_1e_6f = 1e-6f;
}

float max_float(float a, float b) {
	if (a > b)
		return a;
	return b;
}

float min_float(float a, float b) {
	if (a < b)
		return a;
	return b;
}

int max_int(int a, int b) {
	if (a > b)
		return a;
	return b;
}

int min_int(int a, int b) {
	if (a < b)
		return a;
	return b;
}

float my_fabs(float x) {
	if (x < const_float_zero)
		return const_float_minus_one * x;
	return x;
}

float my_sqrt(float x) {
	// TODO implement my own sqrt
	return sqrtf(x);
}

float my_exp(float x) {
	float sum = const_float_one;
	float fact = const_float_one;
	float n = const_float_one;
	char sign = 0;
	if (x < const_float_zero) {
		sign = 1;
		x = const_float_minus_one * x;
	}

	for (int i = 1; i < 50; i++) {
		fact = fact * (x / n);
		sum = sum + fact;
		n = n + const_float_one;
	}

	if (sign)
		return const_float_one / sum;
	return sum;
}

int my_cblas_isamax(const int __N, const float *__X, const int __incX) {
	float vmax = my_fabs(__X[0]);
	int idx = 0;
	for (int i = __incX; i < __N; i += __incX) {
		if (my_fabs(__X[i]) > vmax) {
			vmax = my_fabs(__X[i]);
			idx = i;
		}
	}
	return idx;
}

float my_cblas_sasum(const int __N, const float *__X, const int __incX) {
	float sum = my_fabs(__X[0]);
	for (int i = __incX; i < __N; i += __incX)
		sum = sum + my_fabs(__X[i]);
	return sum;
}

void my_cblas_saxpy(const int __N, const float __alpha, const float *__X, const int __incX, float *__Y, const int __incY) {
	// we suppose the strides are identical
	CHECK_EQ(__incX, __incY);
#pragma omp parallel for shared(__N, __incY, __X, __Y, __alpha) private(i)
	for (int i = 0; i < __N; i += __incY)
		__Y[i] = __Y[i] + __alpha * __X[i];
}

float my_cblas_sdot(const int __N, const float *__X, const int __incX, const float *__Y, const int __incY) {
	// we suppose the strides are identical
	CHECK_EQ(__incX, __incY);
	float sum = __X[0] * __Y[0];
	for (int i = __incY; i < __N; i += __incY) {
		sum = sum + __X[i] * __Y[i];
	}
	return sum;
}

void my_cblas_sscal(const int __N, const float __alpha, float *__X, const int __incX) {
#pragma omp parallel for shared(__N, __incX, __X, __alpha) private(i)
	for (int i = 0; i < __N; i += __incX) {
		__X[i] = __X[i] * __alpha;
	}
}

float my_cblas_snrm2(const int __N, const float *__X, const int __incX) {
	float norm = __X[0] * __X[0];
	for (int i = __incX; i < __N; i += __incX)
		norm = norm + __X[i] * __X[i];
	return my_sqrt(norm);
}

void CxMxM(const float* A, const float* B, float* C, size_t M, size_t N, size_t K, float alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = const_float_zero;
			float sum = const_float_zero;
			for (k = 0; k < K; k++)
				sum += A[i*K+k] * B[k*N+j];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMTxM(const float* A, const float* B, float* C, size_t M, size_t N, size_t K, float alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = const_float_zero;
			float sum = const_float_zero;
			for (k = 0; k < K; k++)
				sum += A[k*M+i] * B[k*N+j];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMxMT(const float* A, const float* B, float* C, size_t M, size_t N, size_t K, float alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = const_float_zero;
			float sum = const_float_zero;
			for (k = 0; k < K; k++)
				sum += A[i*K+k] * B[j*K+k];
			C[i*N+j] += alpha * sum;
		}
	}
}

void CxMTxMT(const float* A, const float* B, float* C, size_t M, size_t N, size_t K, float alpha, char reset) {
	size_t i, j, k;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (reset)
				C[i*N+j] = const_float_zero;
			float sum = const_float_zero;
			for (k = 0; k < K; k++)
				sum += A[k*M+i] * B[j*K+k];
			C[i*N+j] += alpha * sum;
		}
	}
}

void MpM(float* A, const float* B, size_t M, size_t N) {
	size_t i, j;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			A[i*N+j] += B[i*N+j];
		}
	}
}

void CxM(float* A, const float alpha, size_t M, size_t N) {
	size_t i, j;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			A[i*N+j] *= alpha;
}

void my_cblas_sgemv(const enum MY_CBLAS_ORDER __Order,
		const enum MY_CBLAS_TRANSPOSE __TransA,
		const int __M, const int __N,
		const float __alpha,
		const float *__A, const int __lda,
		const float *__X, const int __incX,
		const float __beta,
		float *__Y, const int __incY) {

	// LOG(INFO) << "In my_cblas_sgemv " << (__TransA == CblasTrans);

	int n = (__TransA == MyCblasTrans) ? __M : __N;
	int m = (__TransA == MyCblasTrans) ? __N : __M;

	if (__TransA == MyCblasTrans) {
		for (int r = 0; r < m; r++) {
			float sum = const_float_zero;
			for (int c = 0; c < n; c++) {
				int idx = c * m + r;
				sum += __A[idx] * __X[c];
			}
			__Y[r] = __alpha * sum + __beta * __Y[r];
		}
	}
	else {
		for (int r = 0; r < m; r++) {
			float sum = const_float_zero;
			for (int c = 0; c < n; c++) {
				int idx = r * n + c;
				sum += __A[idx] * __X[c];
			}
			__Y[r] = __alpha * sum + __beta * __Y[r];
		}
	}
}

void my_cblas_sgemm(const enum MY_CBLAS_ORDER __Order,
		const enum MY_CBLAS_TRANSPOSE __TransA,
		const enum MY_CBLAS_TRANSPOSE __TransB,
		const int __M, const int __N, const int __K,
		const float __alpha,
		const float *__A, const int __lda,
		const float *__B, const int __ldb,
		const float __beta,
		float *__C, const int __ldc) {

	// LOG(INFO) << "In my_cblas_sgemm";

	if (__TransA == MyCblasNoTrans) {
		if (__TransB == MyCblasNoTrans) {
			CxM(__C, __beta, __M, __N);
			CxMxM(__A, __B, __C, __M, __N, __K, __alpha, 0);
		}
		else {
			CxM(__C, __beta, __M, __N);
			CxMxMT(__A, __B, __C, __M, __N, __K, __alpha, 0);
		}
	}
	else {
		if (__TransB == MyCblasNoTrans) {
			CxM(__C, __beta, __N, __K);
			CxMTxM(__A, __B, __C, __M, __N, __K, __alpha, 0);
		}
		else {
			CxM(__C, __beta, __N, __N);
			CxMTxMT(__A, __B, __C, __M, __N, __K, __alpha, 0);
		}
	}
}

} // end namespace singa
