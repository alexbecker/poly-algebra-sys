// matrices.c
// used for calculating determinants and inverses

#include "matrices.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#define alt(i) ((i % 2 == 0) ? 1 : -1)
#define EFF_DET_CUTOFF 5

/* ---------- Constructors ---------- */

// returns the n-by-n identity matrix
matrix *identity_matrix(int n) {
	matrix *result = alloc_matrix(n, n);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			if (i == j) {
				result->entries[i][j] = 1.0;
			} else {
				result->entries[i][j] = 0.0;
			}
		}
	}
	
	return result;
}

/* ---------- Memory Functions ---------- */

matrix *alloc_matrix(int m, int n) {
	matrix *result = (matrix*) malloc(sizeof(matrix));
	result->m = m;
	result->n = n;
	result->entries = (vector*) malloc(sizeof(vector) * m);
	for (int i=0; i<m; i++) {
		result->entries[i] = (vector) malloc(sizeof(matrix_entry) * n);
	}
	
	return result;
}

void free_matrix(matrix *a) {
	for (int i=0; i<a->m; i++) {
		free(a->entries[i]);
	}
	free(a->entries);
	free(a);
}

matrix *copy_matrix(matrix a) {
	matrix *cpy = (matrix*) malloc(sizeof(matrix));
	cpy->m = a.m;
	cpy->n = a.n;
	cpy->entries = (vector*) malloc(sizeof(vector) * cpy->m);
	for (int i=0; i<cpy->m; i++) {
		cpy->entries[i] = (vector) malloc(sizeof(matrix_entry) * cpy->n);
		for (int j=0; j<cpy->n; j++) {
			cpy->entries[i][j] = a.entries[i][j];
		}
	}
	
	return cpy;
}

/* ---------- Matrix Arithmetic ---------- */

vector matrix_eval(matrix a, vector v) {
	vector result = (vector) calloc(a.m, sizeof(matrix_entry));
	for (int i=0; i<a.m; i++) {
		for (int j=0; j<a.n; j++) {
			result[i] += a.entries[i][j] * v[j];
		}
	}
	
	return result;
}

// runtime: O(a.m*a.n*b.n)
matrix *matrix_mult(matrix a, matrix b) {
	assert(a.n == b.m);	// otherwise product is not defined
	matrix *result = alloc_matrix(a.m, b.n);
	for (int i=0; i<result->m; i++) {
		for (int j=0; j<result->n; j++) {
			result->entries[i][j] = 0;
			for (int k=0; k<a.n; k++)
				result->entries[i][j] += a.entries[i][k] * b.entries[k][j];
		}
	}
	return result;
}

// returns a minus the zeroth row and nth column
static matrix *matrix_minor(matrix a, int n) {
	matrix *result = alloc_matrix(a.m - 1, a.n - 1);
	for (int i=0; i<result->m; i++) {
		for (int j=0; j<result->n; j++) {
			if (j < n) {
				result->entries[i][j] = a.entries[i+1][j];
			} else {
				result->entries[i][j] = a.entries[i+1][j+1];
			}
		}
	}

	return result;
}

// inefficient determinant algorithm
// runtime: O(n!)
static matrix_entry ineff_det(matrix a) {
	if (a.n == 1)
		return a.entries[0][0];
	matrix_entry result = 0;
	for (int i=0; i<a.n; i++) {
		result += ((matrix_entry) alt(i)) * a.entries[0][i] * ineff_det(*matrix_minor(a, i));
	}
	
	return result;
}

// exchanges row1 and row2 in a
static void exchange_rows(matrix a, int row1, int row2) {
	printf("exchange\n");
	assert((a.m > row1) && (a.m > row2));
	for (int i=0; i<a.n; i++) {
		matrix_entry temp_entry = a.entries[row1][i];
		a.entries[row1][i] = a.entries[row2][i];
		a.entries[row2][i] = temp_entry;
	}
}

// given a matrix a, converts it to a matrix u and returns p and l where pa = lu
// adds number of rows exchanged to row_exchanges (if unneeded just use NULL)
// uses the doolittle algorithm
// runtime: O(n^3)
static matrix **lu_decomp(matrix *a, int *row_exchanges) {
	assert(a->n == a->m);	// only defined for square matrices
	matrix **pl = (matrix**) malloc(sizeof(matrix*) * 2);
	pl[0] = identity_matrix(a->m);	// the matrix p
	pl[1] = identity_matrix(a->m);	// the matrix l

	for (int n=0; n<a->m-1; n++) {
		// permute a so that a->entries[n][n] is invertible
		if (fabs(a->entries[n][n]) < DBL_MIN) {
			// find a row below n with nonzero nth entry
			int new_row = n;
			while (fabs(a->entries[++new_row][n]) < DBL_MIN) {
				// assert that we are not at the last row
				assert(new_row < a->m);
			}
			// exchange appropriate rows of a and p
			exchange_rows(*a, n, new_row);
			exchange_rows(*pl[0], n, new_row);
			if (row_exchanges != NULL)
				*row_exchanges++;
		}
		// perform elimination
		for (int i=n+1; i<a->m; i++) {
			matrix_entry l_entry = a->entries[i][n] / a->entries[n][n];
			pl[1]->entries[i][n] = l_entry; // update l
			// subtract l_entry * the nth row of a from the ith row
			for (int j=0; j<a->n; j++) {
				a->entries[i][j] -= l_entry * a->entries[n][j];
			}
		}
	}
	
//	pl[0] = transpose(pl[0]); // inverts p, giving a = plu
	return pl;
}

// determinant computed using LU-decomposition
// runtime: O(n^3)
static matrix_entry eff_det(matrix a) {
	// make a copy of a so that the entries are not affected by lu_decomp
	matrix *a_cpy = copy_matrix(a);
	
	int row_exchanges = 0;
	lu_decomp(a_cpy, &row_exchanges);
	// a is now upper triangular with same determinant up to sign
	matrix_entry result = alt(row_exchanges);	// takes care of sign
	for (int i=0; i<a_cpy->m; i++) {
		result *= a_cpy->entries[i][i];
	}
	free_matrix(a_cpy);
	 
	return result;
}

matrix_entry det(matrix a) {
	assert(a.m == a.n);	// only defined for square matrices
	// choose which determinant algorithm is faster
	if (a.m < EFF_DET_CUTOFF)
		return ineff_det(a);
	return eff_det(a);
}

static matrix *transpose(matrix a) {
	matrix *result = alloc_matrix(a.n, a.m);
	for (int i=0; i<result->m; i++) {
		for (int j=0; j<result->n; j++)
			result->entries[i][j] = a.entries[j][i];
	}
	return result;
}

// inverts a lower-triangular matrix by forward substitution
// runtime: O(n^3)
static matrix *invert_lower_tri_matrix(matrix a) {
	matrix *identity = identity_matrix(a.m);
	// first we calculate the transpose of the result, as this is easier
	matrix *result_transposed = alloc_matrix(a.m, a.n);
	for (int i=0; i<a.m; i++) {
		for (int j=0; j<a.n; j++) {
			matrix_entry sum_previous_entries = 0;
			for (int k=0; k<j; k++) {
				sum_previous_entries += a.entries[j][k] * result_transposed->entries[i][k];
			}
			result_transposed->entries[i][j] = (identity->entries[i][j] - sum_previous_entries) / a.entries[j][j];
		}
	}
	free_matrix(identity);
	matrix *result = transpose(*result_transposed);
	free_matrix(result_transposed);
	
	return result;
}

// inverts an upper-triangular matrix by converting to a lower-triangular one and using the above algorithm
// runtime: O(n^3)
static matrix *invert_upper_tri_matrix(matrix a) {
	matrix *identity = identity_matrix(a.m);
	matrix *lower_tri = transpose(a);
	// same as before, only we do not transpose the result
	matrix *result = alloc_matrix(a.m, a.n);
	for (int i=0; i<a.m; i++) {
		for (int j=0; j<a.n; j++) {
			matrix_entry sum_previous_entries = 0;
			for (int k=0; k<j; k++) {
				sum_previous_entries += lower_tri->entries[j][k] * result->entries[i][k];
			}
			result->entries[i][j] = (identity->entries[i][j] - sum_previous_entries) / lower_tri->entries[j][j];
		}
	}
	free_matrix(identity);
	free_matrix(lower_tri);
	
	return result;
}

// runtime: O(n^3)
matrix *invert_matrix(matrix a) {
	assert(a.m == a.n);	// only defined for square matrices
	// make a copy of a so as not to modify it
	matrix *u = copy_matrix(a);
	// perform lu decomposition
	matrix **pl = lu_decomp(u, NULL);
	matrix *u_inv = invert_upper_tri_matrix(*u);
	matrix *l_inv = invert_lower_tri_matrix(*pl[1]);
	matrix *lu_inv = matrix_mult(*u_inv, *l_inv);
	matrix *result = matrix_mult(*lu_inv, *pl[0]);
	free_matrix(u);
	free_matrix(pl[0]);
	free_matrix(pl[1]);
	free_matrix(u_inv);
	free_matrix(l_inv);
	free_matrix(lu_inv);
	free(pl);
	
	return result;
}

/* ---------- Input / Output ---------- */

// prints a matrix to stdout
void print_matrix(matrix a) {
	for (int i=0; i<a.m; i++) {
		for (int j=0; j<a.n; j++) {
			printf("%lf ", (double) a.entries[i][j]);
		}
		printf("\n");
	}
}

/* ---------- Testing ---------- */
// to test, run "make test matrices"

#ifdef TEST_MATRICES

/* should output:
 * print_matrix:
 *  3  1 -1
 *  0 -2  0
 *  5  0  0
 * matrix minor:
 *  0  0
 *  5  0
 * ineff_det: -10
 * eff_det: 10 
 * lu decomposition is not unique, output can be checked manually
 * invert_lower_tri_matrix: 
 * 1.000000 0.000000 0.000000 
 * 0.000000 1.000000 0.000000 
 * -1.666667 -0.833333 1.000000 
 * invert_upper_tri_matrix: 
 * 0.333333 0.166667 0.200000 
 * 0.000000 -0.500000 0.000000 
 * 0.000000 -0.000000 0.600000 
 * invert_matrix: 
 * 0.000000 0.000000 0.200000 
 * 0.000000 -0.500000 0.000000 
 * -1.000000 -0.500000 0.600000 */
void test_matrix_functions() {
	matrix *a = identity_matrix(3);
	a->entries[0][0] = 3.0;
	a->entries[0][1] = 1.0;
	a->entries[0][2] = -1.0;
	a->entries[1][1] = -2.0;
	a->entries[2][0] = 5.0;
	a->entries[2][2] = 0.0;
	printf("print_matrix: \n");
	print_matrix(*a);
	printf("matrix_minor: \n");
	print_matrix(*matrix_minor(*a, 1));
	printf("ineff_det: %lf\n", (double) ineff_det(*a));
	printf("eff_det: %lf\n", (double) eff_det(*a));
	matrix *a_cpy = copy_matrix(*a);
	printf("lu_decomp: \n");
	matrix **pl = lu_decomp(a_cpy, (int*) NULL);
	print_matrix(*pl[0]);	// prints p
	print_matrix(*pl[1]);	// prints l
	print_matrix(*a_cpy);	// prints u
	printf("invert_lower_tri_matrix: \n");
	print_matrix(*invert_lower_tri_matrix(*pl[1]));
	printf("invert_upper_tri_matrix: \n");
	print_matrix(*invert_upper_tri_matrix(*a_cpy));
	printf("invert_matrix: \n");
	print_matrix(*invert_matrix(*a));
}

int main(int argc, char **argv) {
	test_matrix_functions();
	exit(0);
}

#endif
