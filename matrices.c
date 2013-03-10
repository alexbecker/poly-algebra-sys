// matrices.c
// used for calculating determinants

#include "matrices.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define alt(i) ((i % 2 == 0) ? 1 : -1)

static matrix *copy_matrix(matrix a) {
	matrix *cpy = (matrix*) malloc(sizeof(matrix));
	cpy->m = a.m;
	cpy->n = a.n;
	cpy->entries = (vector*) malloc(sizeof(vector) * cpy->m);
	for (int i=0; i<cpy->m; i++) {
		cpy->entries[i] = (vector) malloc(sizeof(entry_type) * cpy->n);
		for (int j=0; j<cpy->n; j++) {
			cpy->entries[i][j] = a.entries[i][j];
		}
	}
	
	return cpy;
}

// returns a minus the zeroth row and nth column
static matrix *matrix_minor(matrix a, int n) {
	matrix *result = (matrix*) malloc(sizeof(matrix));
	result->m = a.n - 1;
	result->n = a.m - 1;
	result->entries = (vector*) malloc(sizeof(vector) * result->m);
	for (int i=0; i<result->m; i++) {
		result->entries[i] = (vector) malloc(sizeof(entry_type) * result->n);
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
entry_type ineff_det(matrix a) {
	assert(a.n == a.m);	// only defined for square matrices
	if (a.n == 1)
		return a.entries[0][0];
	entry_type result = 0;
	for (int i=0; i<a.n; i++) {
		result += ((entry_type) alt(i)) * a.entries[0][i] * ineff_det(*matrix_minor(a, i));
	}
	
	return result;
}

// exchanges row1 and row2 in a
static void exchange_rows(matrix a, int row1, int row2) {
	assert((a.m > row1) && (a.m > row2));
	for (int i=0; i<a.n; i++) {
		entry_type temp_entry = a.entries[row1][i];
		a.entries[row1][i] = a.entries[row2][i];
		a.entries[row2][i] = temp_entry;
	}
}

// given a matrix a, converts it to a matrix u and returns p and l where pa = lu
// adds number of rows exchanged to row_exchanges (if unneeded just use NULL)
// uses the doolittle algorithm
// runtime: O(n^3)
static matrix *lu_decomp(matrix a, int *row_exchanges) {
	assert(a.n == a.m);	// only defined for square matrices
	matrix *pl = (matrix*) malloc(sizeof(matrix) * 2);	// pl[0] is p, pl[1] is l
	pl[0].m = a.m;
	pl[0].n = a.n;
	pl[1].m = a.m;
	pl[1].n = a.n;
	pl[0].entries = (vector*) malloc(sizeof(vector) * a.m);
	pl[1].entries = (vector*) malloc(sizeof(vector) * a.m);
	for (int i=0; i<a.m; i++) {
		pl[0].entries[i] = (vector) calloc(a.n, sizeof(entry_type));
		pl[1].entries[i] = (vector) calloc(a.n, sizeof(entry_type));
		pl[0].entries[i][i] = (entry_type) 1;
		pl[1].entries[i][i] = (entry_type) 1;
	}
	// main part of algorithm
	for (int n=0; n<a.m-1; n++) {
		// permute a so that a->entries[n][n] is invertible
		if (a.entries[n][n] == (entry_type) 0) {
			int new_row = n;
			// find a row below n with nonzero nth entry
			while (a.entries[++new_row][n] == (entry_type) 0) {}
			// exchange appropriate rows of a and p
			exchange_rows(a, n, new_row);
			exchange_rows(pl[0], n, new_row);
			if (row_exchanges != NULL)
				*row_exchanges++;
		}
		// perform elimination
		for (int i=n+1; i<a.m; i++) {
			entry_type l_entry = a.entries[i][n] / a.entries[n][n];
			pl[1].entries[i][n] = l_entry; // update l
			// subtract l_entry * the nth row of a from the ith row
			for (int j=0; j<a.n; j++) {
				a.entries[i][j] -= l_entry * a.entries[n][j];
			}
		}
	}
	
//	pl[0] = *transpose(pl[0]); // inverts p, giving a = plu
	return pl;
}

// determinant computed using LU-decomposition
// runtime: O(n^3)
entry_type eff_det(matrix a) {
	assert(a.m == a.n);	// only defined for square matrices]

	// make a copy of a so that the entries are not affected by lu_decomp
	matrix a_cpy = *copy_matrix(a);
	
	int row_exchanges = 0;
	lu_decomp(a_cpy, &row_exchanges);
	// a is no upper triangular with same determinant up to sign
	entry_type result = (entry_type) alt(row_exchanges);	// takes care of sign
	for (int i=0; i<a_cpy.m; i++)
		result *= a_cpy.entries[i][i];
	return result;
}

// prints a matrix to stdout
// used primarily for testing
void print_matrix(matrix a) {
	for (int i=0; i<a.m; i++) {
		for (int j=0; j<a.n; j++) {
			printf("%lf ", (double) a.entries[i][j]);
		}
		printf("\n");
	}
}

// used only for testing
// to test, run "make test matrices"
#ifdef TEST_MATRICES

// should output:
// print_matrix:
//  3  1 -1
//  0 -2  0
//  5  0  0
// matrix minor:
//  0  0
//  5  0
// ineff_det: -10
// eff_det: 10
// lu decomposition output can be checked manually
void test_matrix_functions() {
	matrix *a = (matrix*) malloc(sizeof(matrix));
	a->m = 3;
	a->n = 3;
	a->entries = (vector*) malloc(sizeof(vector) * a->m);
	for (int i=0; i<a->m; i++) {
		a->entries[i] = (vector) calloc(a->n, sizeof(entry_type));
	}
	a->entries[0][0] = (entry_type) 3;
	a->entries[0][1] = (entry_type) 1;
	a->entries[0][2] = (entry_type) -1;
	a->entries[1][1] = (entry_type) -2;
	a->entries[2][0] = (entry_type) 5;
	printf("print_matrix: \n");
	print_matrix(*a);
	printf("matrix_minor: \n");
	print_matrix(*matrix_minor(*a, 1));
	printf("ineff_det: %lf\n", (double) ineff_det(*a));
	printf("eff_det: %lf\n", (double) eff_det(*a));
	printf("lu_decomp: \n");
	matrix *pl = lu_decomp(*a, (int*) NULL);
	print_matrix(pl[0]);	// prints p
	print_matrix(pl[1]);	// prints l
	print_matrix(*a);		// prints u
}

int main(int argc, char **argv) {
	test_matrix_functions();
	exit(0);
}

#endif
