// algebraics.h

#ifndef ALGEBRAICS_H
#define ALGEBRAICS_H

#include "roots.h"
#include "precision.h"

typedef struct algebraic {
	ball approx_val;
	polynomial *minimal_polynomial;	// a NULL pointer indicates that the minimal polynomial is not known
} algebraic;

void free_algebraic(algebraic *);

algebraic *copy_algebraic(algebraic);

void define_minimal_polynomial(algebraic *, int, int);

int is_uniquely_defined(algebraic);

void refine_approx_val(algebraic *, root_type);

algebraic *add_algebraics(algebraic, algebraic);

algebraic *negate_algebraic(algebraic);

algebraic *subtract_algebraics(algebraic, algebraic);

algebraic *mult_algebraics(algebraic, algebraic);

algebraic *invert_algebraic(algebraic);

algebraic *divide_algebraics(algebraic, algebraic);

algebraic *read_algebraic();

void print_algebraic(algebraic);

void print_all_information(algebraic);

#endif
