// matrices.h

#ifndef MATRICES_H
#define MATRICES_H

#include "precision.h"

typedef entry_type* vector;

typedef struct matrix {
	int m, n;	// m rows, n columns
	vector *entries;	// entries[i] is ith row
} matrix;

entry_type ineff_det(matrix a);

entry_type eff_det(matrix a);

void print_matrix(matrix a);

#endif
