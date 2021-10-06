#include <iostream>
#include "singa/core/posit.h"

namespace singa {

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

void compare_posit(uint32_t a, uint32_t b, int ps, int es);
void unpack_exp(int exp, int es, int* r, int* e);

void pack_posit(posit_t p, uint32_t* posit, int* err);
void unpack_posit(int ps, int es, uint32_t posit, int* special, int* sign, int* k, int* e, int* exp, int* fs, int* f, int* err);
uint32_t encode_posit(int ps, int es, double val, int* err);
double decode_posit(int ps, int es, uint32_t posit, int* err);


}