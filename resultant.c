// resultant.c
// uses the resultant to, given polynomials p and q, find polynomials
// with roots at the sums and products of the roots of p and q

#include "resultant.h"
#include "matrices.h"
#include "interpolate.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// returns the sylvester matrix w.r.t. y of p(y) and q(x - y) evaluated at x = x_val
static matrix *sylvester_matrix_sum(polynomial p, polynomial q, int x_val) {
	assert((p.deg > 0) && (q.deg > 0));	// the sylvester matrix is not defined for constant polynomials
	// get the coefficients of q(x_val - y)
	polynomial *q_shifted = linear_change_of_variables(q, -1, x_val);
	// calculate sylvester matrix
	matrix *result = alloc_matrix(p.deg + q_shifted->deg, p.deg + q_shifted->deg);
	for (int i=0; i<q_shifted->deg; i++) {
		for (int j=0; j<i; j++) {
			result->entries[i][j] = 0;
		}
		for (int j=i; j<=i+p.deg; j++) {
			result->entries[i][j] = p.coefficients[j - i];
		}
		for (int j=i+p.deg+1; j<p.deg+q_shifted->deg; j++) {
			result->entries[i][j] = 0;
		}
	}
	for (int i=q_shifted->deg; i<p.deg+q_shifted->deg; i++) {
		for (int j=0; j<i-q_shifted->deg; j++) {
			result->entries[i][j] = 0;
		}
		for (int j=i-q_shifted->deg; j<=i; j++) {
			result->entries[i][j] = q_shifted->coefficients[j - i + q_shifted->deg];
		}
		for (int j=i+1; j<p.deg+q_shifted->deg; j++) {
			result->entries[i][j] = 0;
		}
	}
	// free intermediates
	free_polynomial(q_shifted);
	
	return result;
}

static matrix_entry scalar_resultant_sum(polynomial p, polynomial q, matrix_entry x_val) {
	matrix *sylv_matrix = sylvester_matrix_sum(p, q, x_val);
	matrix_entry result = det(*sylv_matrix);
	free_matrix(sylv_matrix);
	
	return result;
}

// computes a polynomial with zeros at the sums of the zeros of p and q
polynomial *resultant_sum(polynomial p, polynomial q) {
	// evaluate the resultant at 0,...,p.deg*q.deg
	vector vals = (vector) malloc(sizeof(matrix_entry) * (p.deg * q.deg + 1));
	for (int i=0; i<=p.deg*q.deg; i++) {
		vals[i] = scalar_resultant_sum(p, q, i);
	}
	// interpolate to recover the desired polynomial
	polynomial *result = interpolate(vals, p.deg * q.deg);
	free(vals);
	
	return result;
}

// returns the sylvester matrix w.r.t. y of p(y) and y^n*q(x/y) evaluated at x = x_val
static matrix *sylvester_matrix_product(polynomial p, polynomial q, int x_val) {
	assert((p.deg > 0) && (q.deg > 0));	// the sylvester matrix is not defined for constant polynomials
	// get the coefficients of y^n*q(x_val/y)
	polynomial *q_fixed_x = linear_change_of_variables(q, x_val, 0);
	polynomial *q_reversed = reverse_polynomial(*q_fixed_x);
	free_polynomial(q_fixed_x);
	// calculate sylvester matrix
	matrix *result = alloc_matrix(p.deg + q_reversed->deg, p.deg + q_reversed->deg);
	for (int i=0; i<q_reversed->deg; i++) {
		for (int j=0; j<i; j++) {
			result->entries[i][j] = 0;
		}
		for (int j=i; j<=i+p.deg; j++) {
			result->entries[i][j] = p.coefficients[j - i];
		}
		for (int j=i+p.deg+1; j<p.deg+q_reversed->deg; j++) {
			result->entries[i][j] = 0;
		}
	}
	for (int i=q_reversed->deg; i<p.deg+q_reversed->deg; i++) {
		for (int j=0; j<i-q_reversed->deg; j++) {
			result->entries[i][j] = 0;
		}
		for (int j=i-q_reversed->deg; j<=i; j++) {
			result->entries[i][j] = q_reversed->coefficients[j - i + q_reversed->deg];
		}
		for (int j=i+1; j<p.deg+q_reversed->deg; j++) {
			result->entries[i][j] = 0;
		}
	}
	// free intermediates
	free_polynomial(q_reversed);
	
	return result;
}

static matrix_entry scalar_resultant_product(polynomial p, polynomial q, matrix_entry x_val) {
	matrix *sylv_matrix = sylvester_matrix_product(p, q, x_val);
	matrix_entry result = det(*sylv_matrix);
	free_matrix(sylv_matrix);
	
	return result;
}

// computes a polynomial with zeros at the sums of the products of p and q
polynomial *resultant_product(polynomial p, polynomial q) {
	// evaluate the resultant at 0,...,p.deg*q.deg
	vector vals = (vector) malloc(sizeof(matrix_entry) * (p.deg * q.deg + 1));
	// determinant algorithm doesn't always work for determinant 0, so this is kind of a hack
	vals[0] = 0;
	for (int i=1; i<=p.deg*q.deg; i++) {
		vals[i] = scalar_resultant_product(p, q, i);
	}
	// interpolate to recover the desired polynomial
	polynomial *result = interpolate(vals, p.deg * q.deg);
	free(vals);
	
	return result;
}

/* ---------- Testing ---------- */
// to test, run "make test resultant"

#ifdef TEST_RESULTANT

/* should output:
 * resultant_sum: -8x^6 - 16x^5 + 112x^4 + 256x^3 - 168x^2 - 496x^1 - 192
 * x^10 - 10x^8 + 38x^6 + 2x^5 - 100x^4 + 40x^3 + 121x^2 + 38x^1 - 17
 * resultant_product: 8x^6 + 32x^5 - 160x^4 + 128x^3 */
void test_resultant_functions() {
	polynomial *p = alloc_polynomial(3);
	p->coefficients[0] = -1;
	p->coefficients[1] = 2;
	p->coefficients[2] = 5;
	p->coefficients[3] = 2;
	polynomial *q = alloc_polynomial(2);
	q->coefficients[0] = -2;
	q->coefficients[1] = -4;
	q->coefficients[2] = 0;
	polynomial *f = calloc_polynomial(5);
	f->coefficients[0] = 1;
	f->coefficients[4] = -1;
	f->coefficients[5] = 1;
	polynomial *g = calloc_polynomial(2);
	g->coefficients[0] = 1;
	g->coefficients[2] = -2;
	printf("resultant_sum: ");
	print_polynomial(*resultant_sum(*p, *q));
	print_polynomial(*resultant_sum(*f, *g));
	printf("resultant_product: ");
	print_polynomial(*resultant_product(*p, *q));
}

int main(int argc, char **argv) {
	test_resultant_functions();
	exit(0);
}

#endif
