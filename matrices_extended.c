// graveyard for code excluded from matrices.c

// runtime: O(a.m*a.n*b.n)
matrix *matrix_mult(matrix a, matrix b) {
	assert(a.n == b.m);	// otherwise product is not defined
	matrix *result = (matrix*) malloc(sizeof(matrix));
	result.m = a.m;
	result.n = b.n;
	result.entries = (vector*) malloc(sizeof(vector) * result.m);
	for (int i=0; i<result.m; i++) {
		result.entries[i] = (vector) malloc(sizeof(entry_type) * result.n);
		for (int j=0; j<result.n; j++) {
			result.entries[i][j] = 0;
			for (int k=0; k<a.n; k++)
				result.entries[i][j] += a.entries[i][k] * b.entries[k][j];
		}
	}
	return result;
}

double vec_norm(vector v, int length) {
	entry_type result = 0;
	for (int i=0; i<length; i++)
		result += v[i] * v[i];
	return sqrt((double) result);
}

static matrix *transpose(matrix a) {
	matrix *result = (matrix*) malloc(sizeof(matrix));
	result->m = a.n;
	result->n = a.m;
	result->entries = (vector*) malloc(sizeof(vector) * result->m);
	for (int i=0; i<result->m, i++) {
		result->entries[i] = (vector) malloc(sizeof(entry_type) * result->n);
		for (int j=0; j<result->n; j++)
			result->entries[i][j] = a.entries[j][i];
	}
	return result;
}

// returns a random matrix with double entries between -1 and 1
// runtime O(mn)
matrix *rand_matrix(int m, int n) {
	matrix *result = (matrix*) malloc(sizeof(matrix));
	result->m = m;
	result->n = n;
	srand(time(NULL)); // initialize random seed
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++)
			result->entries[i][j] = ((entry_type) alt(rand())) * ((entry_type) rand()) / ((entry_type) RAND_MAX);
	}
	return result;
}

// UNFINISHED
matrix *matrix_inv(matrix a) {
	matrix *result = (matrix) malloc(sizeof(matrix));

// REQUIRES matrix_inv
// given a matrix a, returns a matrix which represents orthogonal projection onto the columns of a
matrix *orth_proj_matrix(matrix a) {
	matrix a_trans = *transpose(a);
	return matrix_mult(a, *matrix_mult(*matrix_inv(*matrix_mult(a_trans, a)), a_trans));
}
