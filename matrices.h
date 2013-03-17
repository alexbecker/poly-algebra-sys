// matrices.h

#ifndef MATRICES_H
#define MATRICES_H

#include "precision.h"

typedef matrix_entry* vector;

typedef struct matrix {
	int m, n;	// m rows, n columns
	vector *entries;	// entries[i] is ith row
} matrix;

matrix *alloc_matrix(int, int);

void free_matrix(matrix *);

matrix *copy_matrix(matrix);

matrix *identity_matrix(int);

vector matrix_eval(matrix, vector);

matrix_entry det(matrix);

matrix *invert_matrix(matrix);

void print_matrix(matrix);

#endif
