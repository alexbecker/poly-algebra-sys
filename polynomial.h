// polynomial.h

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "precision.h"

// specifically an integer polynomial
typedef struct polynomial {
	int deg;			// degree of zero polynomial is considered 0
	int *coefficients;	// from highest to lowest degree, for algorithmic reasons
						// coefficients[0] MUST be nonzero
} polynomial;

polynomial *int_to_polynomial(int);

polynomial *alloc_polynomial(int);

polynomial *calloc_polynomial(int);

polynomial *copy_polynomial(polynomial);

void free_polynomial(polynomial*);

void strip_leading_zeros(polynomial*);

root_type eval_polynomial(polynomial, root_type);

complex eval_polynomial_complex(polynomial, complex);

polynomial *add_polynomials(polynomial, polynomial);

polynomial *negate_polynomial(polynomial);

polynomial *subtract_polynomials(polynomial, polynomial);

polynomial *increase_degree(polynomial, int);

polynomial *reverse_polynomial(polynomial);

polynomial *scalar_mult_polynomial(polynomial, int);

polynomial *mult_polynomials(polynomial, polynomial);

polynomial *polynomial_power(polynomial, int);

polynomial *compose_polynomials(polynomial, polynomial);

polynomial *linear_change_of_variables(polynomial, int, int);

polynomial *differentiate(polynomial);

polynomial *polynomial_mod(polynomial, polynomial);

polynomial *read_polynomial();

void print_polynomial(polynomial);

#endif
