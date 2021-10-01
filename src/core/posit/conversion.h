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
#ifndef CONVERSION_H_
#define CONVERSION_H_

// define PDEBUG to enable printing info
// #define PDEBUG 1

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

char* get_error(int err);

double decode_posit(int n, int es, uint32_t posit, int* err);

uint32_t encode_posit(int n, int es, double val, int* err);

void pack_posit(int n, int es, int special, int sign, int k, int e, int fs, int f, uint32_t* posit, int* err);

uint32_t big_to_small_posit(int n1, int es1, uint32_t val1, int n2, int es2, int* err);

uint32_t small_to_big_posit(int n1, int es1, uint32_t val1, int n2, int es2, int* err);

#endif /* CONVERSION_H_ */

