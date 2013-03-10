#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "precision.h"

typedef struct polynomial {
	int deg;
	int *coefficients;	// from highest to lowest degree, for algorithmic reasons
} polynomial;

// accounts for error by treating entry_type as a small ball
typedef struct ball {
	entry_type center, radius;
} ball;

entry_type eval_polynomial(polynomial, entry_type);

polynomial *differentiate(polynomial);

// returns 1 if a root exists, 0 if unsure
int test_for_root(polynomial, ball);

// returns an upper bound on the magnitude of the roots of a polynomial
entry_type root_upper_bound(polynomial);

void print_polynomial(polynomial);

#endif
