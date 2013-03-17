// interpolate.c
// contains functions for polynomial interpolation

#include "interpolate.h"
#include "integers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// returns the matrix used to interpolate evaluation at 0,...,degree
// runtime: O(degree^3)
static matrix *interpolation_matrix(int degree) {
	matrix *result_inv = alloc_matrix(degree + 1, degree + 1);
	for (int i=0; i<=degree; i++) {
		for (int j=0; j<=degree; j++) {
			result_inv->entries[i][j] = pow(i, j);
		}
	}
	
	matrix *result = invert_matrix(*result_inv);
	free_matrix(result_inv);
	
	return result;
}

// returns the polynomial which passes through vals at 0,...,degree
// runtime: O(degree^4)
polynomial *interpolate(vector vals, int degree) {
	polynomial *result = alloc_polynomial(degree);
	matrix *interp_matrix = interpolation_matrix(degree);
	vector fractional_coeffs = matrix_eval(*interp_matrix, vals);
	// clear the denominators in fractional_coeffs and reverse entries
	int total_denom = denom(fractional_coeffs[0]);
	for (int i=1; i<=degree; i++) {
		total_denom = lcm(total_denom, denom(fractional_coeffs[i]));
	}
	for (int i=0; i<=degree; i++) {
		result->coefficients[i] = (int) round(total_denom * fractional_coeffs[degree - i]);
	}
	// free intermediates
	free_matrix(interp_matrix);
	free(fractional_coeffs);
	// remove leading zeros from result
	strip_leading_zeros(result);
	
	return result;
}

/* ---------- Testing ---------- */
// to test, run "make test interpolate"

#ifdef TEST_INTERPOLATE

// should output whatever (x-1)(x-2)(x-3)(x-4) is
void test_interpolate() {
	vector vals = (vector) calloc(5, sizeof(matrix_entry));
	vals[0] = 24;
	printf("interpolate: ");
	print_polynomial(*interpolate(vals, 4));
}

int main(int argc, char **argv) {
	test_interpolate();
	exit(0);
}

#endif
