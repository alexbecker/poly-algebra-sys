// polynomial.c
// impliments basic functions for handling integer polynomials

#include "polynomial.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define abs(a) ((a >= 0) ? a : (0 - a))	// only use for integers

// given a polynomial p, makes a copy of p on the heap
static polynomial *copy_polynomial(polynomial p) {
	polynomial *p_cpy = (polynomial*) malloc(sizeof(polynomial));
	p_cpy->deg = p.deg;
	p_cpy->coefficients = (int*) malloc(sizeof(int) * (p.deg + 1));
	memcpy(p_cpy->coefficients, p.coefficients, sizeof(int) * (p.deg + 1));
	
	return p_cpy;
}

// given polynomials p and q, returns the remainder of p mod q
polynomial *polynomial_mod(polynomial p, polynomial q) {
	// check if there is anything to do
	if ((p.deg < q.deg) || (p.deg == 0 && p.coefficients[0] = 0))
		return copy_polynomial(p);	// return a copy of p so as not to
	
	
	polynomial *result = (polynomial*) malloc(sizeof(polynomial));

// evaluate a polynomial using horner's scheme
entry_type eval_polynomial(polynomial p, entry_type x) {
	if (p.deg == 0)
		return (entry_type) p.coefficients[0];
	p.deg--;	// effectively decrease degree of p by 1
	
	return ((entry_type) p.coefficients[p.deg + 1]) + x * eval_polynomial(p, x);
}

// given a polynomial p, returns p'
polynomial *differentiate(polynomial p) {
	polynomial *p_prime = (polynomial*) malloc(sizeof(polynomial));
	p_prime->deg = p.deg - 1;
	p_prime->coefficients = (int*) malloc(sizeof(int) * p.deg);
	for (int i=0; i<p.deg; i++) {
		p_prime->coefficients[i] = (p.deg - i) * p.coefficients[i];
	}
	
	return p_prime;
}

// returns 1 if a root exists, 0 if unsure
int test_for_root(polynomial p, ball b) {
	entry_type left_val = eval_polynomial(p, b.center - b.radius), right_val = eval_polynomial(p, b.center + b.radius);
	if (((left_val < 0 - EVAL_ERR) && (right_val > EVAL_ERR)) || ((left_val > EVAL_ERR) && (right_val < 0 - EVAL_ERR)))
		return 1;
	return 0;
}

// returns an upper bound on the roots of p
// uses multiple bounds and takes the largest
entry_type root_upper_bound(polynomial p) {
	entry_type Fujiwara_bnd, Cauchy_bnd;
	
	Fujiwara_bnd = (entry_type) 2.0 * pow(fabs(((double) p.coefficients[p.deg]) / (2.0 * (double) p.coefficients[0])), 1.0 / ((double) p.deg));
	for (int i=1; i<p.deg; i++) {
		if (Fujiwara_bnd < (entry_type) 2.0 * pow(fabs(((double) p.coefficients[p.deg-i]) / ((double) p.coefficients[0])), 1.0 / ((double) p.deg - i)))
			Fujiwara_bnd = (entry_type) 2.0 * pow(fabs(((double) p.coefficients[p.deg-i]) / ((double) p.coefficients[0])), 1.0 / ((double) p.deg - i));
	}
	
	Cauchy_bnd = (entry_type) 0;
	for (int i=0; i<p.deg; i++) {
		if (Cauchy_bnd < (entry_type) 1.0 + fabs(((double) p.coefficients[p.deg-i]) / ((double) p.coefficients[0])))
			Cauchy_bnd = (entry_type) 1.0 + fabs(((double) p.coefficients[p.deg-i]) / ((double) p.coefficients[0]));
	}
	
	if (Fujiwara_bnd < Cauchy_bnd) {
		return Fujiwara_bnd;
	} else {
		return Cauchy_bnd;
	}
}

// counts the number of roots of a polynomial using Sturm's theorem
// requires p be square-free
int root_cnt(polynomial p) {
	if (p.deg == 0) {
		assert(p.coefficients[0] != 0);
		return 0;
	}
	polynomial p0 = *copy_polynomial(p);
	polynomial p1 = *differentiate(p);
	for (int 

// given a polynomial p and ptr to an int num_roots
// finds the roots of a polynomial using Newton's method
// returns the (ordered) list of real roots and stores the number in num_roots
ball *get_real_roots(polynomial p, int *num_roots) {
	entry_type *roots = (entry_type*) malloc(sizeof(entry_type) * p.deg);

// prints a polynomial to stdout in human-readable form
void print_polynomial(polynomial p) {
	if (p.coefficients[0] == 1) {
		printf("x^%d", p.deg);
	} else if (p.coefficients[0] == -1) {
		printf("-x^%d", p.deg);
	} else if (p.coefficients[0] != 0) {
		printf("%dx^%d", p.coefficients[0], p.deg);
	} 
	for (int i=1; i<p.deg; i++) {
		if (p.coefficients[i] > 0) {
			printf(" + ");
		} else if (p.coefficients[i] < 0) {
			printf(" - ");
		} else {
			continue;
		}
		if (abs(p.coefficients[i]) == 1) {
			printf("x^%d", p.deg - i);
		} else {
			printf("%dx^%d", abs(p.coefficients[i]), p.deg - i);
		}
	}
	if (p.coefficients[p.deg] > 0) {
		printf(" + %d\n", p.coefficients[p.deg]);
	} else if (p.coefficients[p.deg] < 0) {
		printf(" - %d\n", abs(p.coefficients[p.deg]));
	} else {
		printf("\n");
	}
}

// used only for testing
// to test, run "make test polynomial"
#ifdef TEST_POLYNOMIAL
#define ROOT_2 1.414213562373095
#define ROOT_2_ERR 1e-15

// should output:
// print_polynomial: x^2 - 2
// eval_polynomial: 0
// differentiate: 2x^1
// test_for_root: 1
// root_upper_bound: 2
void test_polynomial_functions() {
	polynomial *p = (polynomial*) malloc(sizeof(polynomial));
	p->deg = 2;
	p->coefficients = (int*) malloc(sizeof(int) * 3);
	p->coefficients[0] = 1;
	p->coefficients[1] = 0;
	p->coefficients[2] = -2;
	printf("print_polynomial: ");
	print_polynomial(*p);
	printf("eval_polynomial: %lf\n", (double) eval_polynomial(*p, ROOT_2));
	printf("differentiate: ");
	print_polynomial(*differentiate(*p));
	ball *b = (ball*) malloc(sizeof(ball));
	b->center = ROOT_2;
	b->radius = ROOT_2_ERR;
	printf("test_for_root: %d\n", test_for_root(*p, *b));
	printf("root_upper_bound: %lf\n", (double) root_upper_bound(*p));
}

int main(int argc, char **argv) {
	test_polynomial_functions();
	exit(0);
}

#endif
