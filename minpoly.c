// minpoly.c
// finds the minimal polynomial of a low-degree positive real algebraic number

#include "minpoly.h"
#include "subset_sum.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

// computes integer powers of 2
static int pow2(int n) {
	if (n == 0) {
		return 1;
	} else {
		int half_pow = pow2(n / 2);
		if (n % 2 == 0)
			return half_pow * half_pow;
		return 2 * half_pow * half_pow;
	}
}

// Reduction of Minimal Polynomial to Subset Sum
// given a value x and integers k and deg
// returns a subset sum problem of size k * (deg + 1) 
// if e_i is 0 if the ith element is not part of the certifying subset and 1 if it is, then the polynomial is:
// (1 + e_0*2^0 + ... + e_{k-1}*2^{k-1})
//  + (-2^{k-1} + e_{k}2^0 + ... + e_{2k-1}2^{k-1})x
//  + ...
//  + (-2^{k-1} + e_{k*deg}2^0 + ... + e_{k*(deg+1)-1}2^{k-1})x^deg
static subset_sum_problem *to_subset_sum(entry_type x, int k, int deg) {
	subset_sum_problem *problem = (subset_sum_problem*) malloc(sizeof(subset_sum_problem));
	problem->size = k * (deg + 1);
	problem->items = (entry_type*) malloc(sizeof(entry_type) * problem->size);
	assert(problem->items != NULL);
	// the sum is the negation of the terms in the polynomial which do not depend on e_i
	// this is -1 + 2^(k-1)x + ... + 2^(k-1)x^deg
	// NOTE: breaks if x = 1, but this case is trivial
	problem->sum = (entry_type) -1;
	for (int i=1; i<=deg; i++) {
		problem->sum += ((entry_type) pow2(k-1)) * ((entry_type) pow((double) x, (double) i));
	}
	// items
	for (int i=0; i<=deg; i++) {
		for (int j=0; j<k; j++) {
			problem->items[k * i + j] = (entry_type) ((double) pow2(j) * pow((double) x, (double) i)); 
		}
	}
	return problem;
}

// given a subset output by running subset_sum_certificate on a problem created by to_subset_sum
// returns the associated polynomial
static polynomial *subset_to_polynomial(subset_sum_problem problem, int k, int *subset) {
	polynomial *result = (polynomial*) malloc(sizeof(polynomial));
	result->deg = problem.size / k - 1;
	result->coefficients = (int*) malloc(sizeof(int) * (result->deg + 1));
	
	// determine the coefficients of the result
	// recall that the coefficient of degree n is result->coefficients[result->deg - n]
	for (int i=1; i<=result->deg; i++) {
		result->coefficients[result->deg - i] = 0 - pow2(k - 1);
		for (int j=0; j<k; j++) {
			result->coefficients[result->deg - i] += subset[k*i+j] * pow2(j);
		}
	}
	// constant term
	result->coefficients[result->deg] = 1;
	for (int j=0; j<k; j++) {
		result->coefficients[result->deg] += subset[j] * pow2(j);
	}
	
	return result;
}

// given a ball b and ints max_k and max_deg
// searches for a minimal polynomial with coefficients bounded by 2^{max_k-1} for an element the ball b by increasing k to max_k and degree to max_deg
// returns the polynomial found or NULL on failure
polynomial *find_min_poly(ball b, int max_k, int max_deg) {
	entry_type error;
	polynomial *min_poly;
	for (int i=1; i<=max_deg; i++) {
		subset_sum_problem problem = *to_subset_sum(b.center, max_k, i);
		int *subset = subset_sum_certificate(problem, &error);
		min_poly = subset_to_polynomial(problem, max_k, subset);
		// test if b.center is a root of min_poly
		if (test_for_root(*min_poly, b))
			return min_poly;
	}
	return NULL;
}

// used only for testing
// to test, run "make test minpoly"
#ifdef TEST_MINPOLY
#define ROOT_2 1.414213562373095
#define ROOT_2_ERR 1e-15
#define DEG5_ALG 1.1673039782614186843

// should output a sum of 7+4sqrt(2)
// then 1, 2, 4, sqrt(2), 2sqrt(2), 4sqrt(2), 2, 4, 8
void test_to_subset_sum() {
	printf("to_subset_sum: ");
	subset_sum_problem *problem = to_subset_sum((entry_type) ROOT_2, 3, 2);
	printf("sum: %lf\nitems: ", problem->sum);
	for (int i=0; i<problem->size - 1; i++)
		printf("%lf, ", problem->items[i]);
	printf("%lf\n", problem->items[problem->size-1]);
}

// should output -4x^2 + 8
void test_subset_to_polynomial() {
	subset_sum_problem *problem = to_subset_sum((entry_type) ROOT_2, 3, 2);
	int *subset = (int*) calloc(problem->size, sizeof(int));
	subset[0] = 1;
	subset[1] = 1;
	subset[2] = 1;
	subset[5] = 1;
	printf("subset_to_polynomial: ");
	print_polynomial(*subset_to_polynomial(*problem, 3, subset));
}

// should output -x^2 + 2 and -x^5 + x + 1
void test_find_min_poly() {
	ball *root_2_ball = (ball*) malloc(sizeof(ball));
	root_2_ball->center = ROOT_2;
	root_2_ball->radius = ROOT_2_ERR;
	printf("min_poly for sqrt(2): ");
	print_polynomial(*find_min_poly(*root_2_ball, 3, 2));
	ball *deg5_alg_ball = (ball*) malloc(sizeof(ball));
	deg5_alg_ball->center = DEG5_ALG;
	deg5_alg_ball->radius = 1e-10;
	printf("min_poly for degree 5 algebraic: ");
	print_polynomial(*find_min_poly(*deg5_alg_ball, 3, 5));
}

int main(int argc, char **argv) {
	test_to_subset_sum();
	test_subset_to_polynomial();
	test_find_min_poly();
	exit(0);
}

#endif
