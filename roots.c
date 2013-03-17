// roots.c
// impliments root finding and testing functions for integer polynomials

#include "roots.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

/* ---------- Constructors ---------- */

ball *new_ball(root_type center, root_type radius) {
	ball *result = (ball*) malloc(sizeof(ball));
	result->center = center;
	result->radius = radius;
	
	return result;
}

complex_ball *new_complex_ball(complex center, root_type radius) {
	complex_ball *result = (complex_ball*) malloc(sizeof(complex_ball));
	result->center = center;
	result->radius = radius;
	
	return result;
}

// returns a root_list containing multiplicity copies of the given root
root_list *root_to_root_list(root_type root, root_type error, int multiplicity) {
	root_list *result = alloc_root_list(multiplicity);
	for (int i=0; i<multiplicity; i++) {
		result->roots[i] = *new_ball(root, error);
	}
	
	return result;
}

// preferred to alloc_root_list(0) which may be implimentation-dependent
root_list *empty_root_list() {
	root_list *result = (root_list*) malloc(sizeof(root_list));
	result->num_roots = 0;
	result->roots = (ball*) NULL;
	
	return result;
}

/* ---------- Memory Functions ---------- */

root_list *alloc_root_list(int num_roots) {
	root_list *result = (root_list*) malloc(sizeof(root_list));
	result->num_roots = num_roots;
	result->roots = (ball*) malloc(sizeof(ball) * num_roots);
	
	return result;
}

complex_root_list *alloc_complex_root_list(int num_roots) {
	complex_root_list *result = (complex_root_list*) malloc(sizeof(complex_root_list));
	result->num_roots = num_roots;
	result->roots = (complex_ball*) malloc(sizeof(complex_ball) * num_roots);
	
	return result;
}

ball *copy_ball(ball b) {
	ball *result = (ball*) malloc(sizeof(ball));
	result->center = b.center;
	result->radius = b.radius;
	
	return result;
}

void free_root_list(root_list *r) {
	free(r->roots);
	free(r);
}

void free_complex_root_list(complex_root_list *r) {
	free(r->roots);
	free(r);
}

/* ---------- (Complex) Ball Functions ---------- */

// returns 1 if the given value is in the given ball, 0 otherwise
int test_ball_membership(root_type x, ball b) {
	if ((x > b.center - b.radius) && (x < b.center + b.radius))
		return 1;
	return 0;
}

int test_complex_ball_membership(complex x, complex_ball b) {
	if (cabsl(x - b.center) < b.radius)
		return 1;
	return 0;
}

/* ---------- Root_list Functions ---------- */

// merges two lists of roots
// does not check for repeat elements
root_list *merge_root_lists(root_list a, root_list b) {
	root_list *result = alloc_root_list(a.num_roots + b.num_roots);
	for (int i=0; i<a.num_roots; i++) {
		result->roots[i] = a.roots[i];
	}
	for (int i=0; i<b.num_roots; i++) {
		result->roots[i + a.num_roots] = b.roots[i];
	}
	
	return result;
}

/* ---------- Root Functions ---------- */

// returns 1 if a root exists, 0 if unsure
int test_for_root(polynomial p, ball b) {
	root_type left_val = eval_polynomial(p, b.center - b.radius), right_val = eval_polynomial(p, b.center + b.radius);
	if (((left_val < -EVAL_ERR) && (right_val > EVAL_ERR)) || ((left_val > EVAL_ERR) && (right_val < -EVAL_ERR)))
		return 1;
	return 0;
}

// returns an upper bound on the abs. val. of the roots of p
// uses multiple bounds and takes the smallest
static root_type root_upper_bound(polynomial p) {
	root_type Fujiwara_bnd, Cauchy_bnd;
	
	Fujiwara_bnd =  2.0 * pow(fabs(p.coefficients[p.deg] / (2.0 * (double) p.coefficients[0])), 1.0 / ((double) p.deg));
	for (int i=1; i<p.deg; i++) {
		if (Fujiwara_bnd < 2.0 * pow(fabs(p.coefficients[p.deg-i] / ((double) p.coefficients[0])), 1.0 / ((double) p.deg - i)))
			Fujiwara_bnd = 2.0 * pow(fabs(p.coefficients[p.deg-i] / ((double) p.coefficients[0])), 1.0 / ((double) p.deg - i));
	}
	
	Cauchy_bnd = 0.0;
	for (int i=0; i<p.deg; i++) {
		if (Cauchy_bnd < 1.0 + fabs( p.coefficients[p.deg-i] / ((double) p.coefficients[0])))
			Cauchy_bnd = 1.0 + fabs(p.coefficients[p.deg-i] / ((double) p.coefficients[0]));
	}
	
	if (Fujiwara_bnd < Cauchy_bnd) {
		return Fujiwara_bnd + MIN_ROOT_ERR;
	} else {
		return Cauchy_bnd + MIN_ROOT_ERR;
	}
}

// helper function for root_cnt
// returns negative (p mod q)
static polynomial *negate_polynomial_mod(polynomial p, polynomial q) {
	polynomial *poly_mod = polynomial_mod(p, q);
	polynomial *result = negate_polynomial(*poly_mod);
	free_polynomial(poly_mod);
	
	return result;
}

// counts the number of roots of a polynomial between the bounds using Sturm's theorem
// requires p be square-free
int root_cnt(polynomial p, root_type lower_bnd, root_type upper_bnd) {
	if (p.deg == 0) {
		assert(p.coefficients[0] != 0);
		return 0;
	}
	// create sturm sequence and count sign changes
	polynomial *sturm_seq = (polynomial*) malloc(sizeof(polynomial) * 2);
	sturm_seq[0] = *copy_polynomial(p);
	sturm_seq[1] = *differentiate(p);
	int seq_ind = 2, sgn_changes_lower = 0, sgn_changes_upper = 0, cur_sgn_lower, cur_sgn_upper;
	// initialize the signs of the chain
	// may break at roots of p
	if (eval_polynomial(sturm_seq[0], lower_bnd) > 0) {
		cur_sgn_lower = 1;
	} else {
		cur_sgn_lower = -1;
	}
	if (eval_polynomial(sturm_seq[0], upper_bnd) > 0) {
		cur_sgn_upper = 1;
	} else {
		cur_sgn_upper = -1;
	}
	// take into account the first term of the chain
	root_type lower_val = eval_polynomial(sturm_seq[1], lower_bnd);
	root_type upper_val = eval_polynomial(sturm_seq[1], upper_bnd);
	if ((lower_val > 0) && (cur_sgn_lower == -1)) {
		sgn_changes_lower++;
		cur_sgn_lower = 1;
	} else if ((lower_val < 0) && (cur_sgn_lower == 1)) {
		sgn_changes_lower++;
		cur_sgn_lower = -1;
	}
	if ((upper_val > 0) && (cur_sgn_upper == -1)) {
		sgn_changes_upper++;
		cur_sgn_upper = 1;
	} else if ((upper_val < 0) && (cur_sgn_upper == 1)) {
		sgn_changes_upper++;
		cur_sgn_upper = -1;
	}
	// remaining terms
	while((sturm_seq[0].deg != 0) && (sturm_seq[1].deg != 0)) {
		int *old_coeff_ptr = sturm_seq[seq_ind % 2].coefficients;
		sturm_seq[seq_ind % 2] = *negate_polynomial_mod(sturm_seq[seq_ind % 2], sturm_seq[(seq_ind + 1) % 2]);
		free(old_coeff_ptr);
		
		lower_val = eval_polynomial(sturm_seq[seq_ind % 2], lower_bnd);
		upper_val = eval_polynomial(sturm_seq[seq_ind % 2], upper_bnd);
		if ((lower_val > 0) && (cur_sgn_lower == -1)) {
			sgn_changes_lower++;
			cur_sgn_lower = 1;
		} else if ((lower_val < 0) && (cur_sgn_lower == 1)) {
			sgn_changes_lower++;
			cur_sgn_lower = -1;
		}
		if ((upper_val > 0) && (cur_sgn_upper == -1)) {
			sgn_changes_upper++;
			cur_sgn_upper = 1;
		} else if ((upper_val < 0) && (cur_sgn_upper == 1)) {
			sgn_changes_upper++;
			cur_sgn_upper = -1;
		}
		
		seq_ind++;
	}
	free(sturm_seq[0].coefficients);
	free(sturm_seq[1].coefficients);
	free(sturm_seq);
	return sgn_changes_lower - sgn_changes_upper;
}

// counts all roots using root_cnt
// again requires p be square-free
int total_root_cnt(polynomial p) {
	root_type root_abs_bnd = root_upper_bound(p);
	return root_cnt(p, -root_abs_bnd, root_abs_bnd);
}

// given a polynomial p, upper and lower bounds and an error bound
// returns the list of roots between the bounds with ball radius at most the error
root_list *get_real_roots(polynomial p, root_type lower_bnd, root_type upper_bnd, root_type error) {
	int num_roots = root_cnt(p, lower_bnd, upper_bnd);
	// if there are no roots in the interval, return an empty root_list
	if (num_roots == 0)
		return empty_root_list();
	// if the difference between the average of the bounds and either bound is less than error, we've found a root
	root_type avg = (upper_bnd + lower_bnd) / 2.0, radius = (upper_bnd - lower_bnd) / 2.0;
	if (radius < error)
		return root_to_root_list(avg, radius, num_roots); 
	// else we split the interval in two and combine the results
	root_list *upper_roots = get_real_roots(p, avg, upper_bnd, error);
	root_list *lower_roots = get_real_roots(p, lower_bnd, avg, error);
	root_list *all_roots = merge_root_lists(*lower_roots, *upper_roots);
	free_root_list(lower_roots);
	free_root_list(upper_roots);
	
	return all_roots;
}

root_list *get_all_real_roots(polynomial p, root_type error) {
	root_type root_abs_bnd = root_upper_bound(p);
	return get_real_roots(p, -root_abs_bnd, root_abs_bnd, error);
}

// finds all complex roots via the Durand–Kerner method
// only works if p is square-free
complex_root_list *get_all_roots(polynomial p, root_type error) {
	complex_root_list *result = alloc_complex_root_list(p.deg);
	// initialize the guesses
	complex *root_guesses = (complex*) malloc(sizeof(complex) * p.deg);
	root_guesses[0] = 0.4 * pow(sqrt(p.deg), 1 / ((double) p.deg)) + 0.9 * pow(sqrt(p.deg), 1 / ((double) p.deg)) * I;
	for (int i=1; i<p.deg; i++) {
		root_guesses[i] = root_guesses[0] * root_guesses[i - 1];
	}
	// iteratively improve the guesses until the difference between two iterations is less than half of error
	root_type difference;
	do {
		for (int i=0; i<p.deg; i++) {
			complex newton_denom = 1;
			for (int j=0; j<p.deg; j++) {
				if (i != j)
					newton_denom *= root_guesses[i] - root_guesses[j];
			}
			complex new_root = root_guesses[i] - eval_polynomial_complex(p, root_guesses[i]) / (p.coefficients[0] * newton_denom);
			difference = cabsl(new_root - root_guesses[i]);
			root_guesses[i] = new_root;
		}
	} while (difference >= error / 2);
	// copy the root_guesses to result, along with the error
	for (int i=0; i<p.deg; i++) {
		result->roots[i] = *new_complex_ball(root_guesses[i], error);
	}
	// free intermediates
	free(root_guesses);
	
	return result;
}

/* ---------- Input / Output  ---------- */

// prints a ball to stdout
void print_ball(ball b) {
	printf("Approximate value: %lf, Error: %lf\n", (double) b.center, (double) b.radius);
}

void print_complex_ball(complex_ball b) {
	printf("Approximate value: %lf + %lfi, Error: %lf\n", (double) __real__ b.center, (double) __imag__ b.center, (double) b.radius);
}

// prints a root_list to stdout
void print_root_list(root_list r) {
	for (int i=0; i<r.num_roots; i++) {
		printf("Root %d: ", i);
		print_ball(r.roots[i]);
	}
}

void print_complex_root_list(complex_root_list r) {
	for (int i=0; i<r.num_roots; i++) {
		printf("Root %d: ", i);
		print_complex_ball(r.roots[i]);
	}
}

/* ---------- Testing ---------- */
// to test, run "make test roots"

#ifdef TEST_ROOTS
#define ROOT_2 1.414213562373095
#define ROOT_2_ERR 1e-15
#define DEG5_ERR 1e-10

/* should output:
 * test_for_root: 1
 * negate_polynomial_mod: -4
 * total_root_cnt: 2
 * get_all_real_roots:
 * Approximate value: -1.414214, Error: 0.000000
 * Approximate value: 1.414214, Error: 0.000000 */
void test_root_functions() {
	polynomial p = *alloc_polynomial(2);
	p.coefficients[0] = -1;
	p.coefficients[1] = 0;
	p.coefficients[2] = 2;
	ball root_2_ball = *new_ball(ROOT_2, ROOT_2_ERR);
	printf("test_for_root: %d\n", test_for_root(p, root_2_ball));
	printf("negate_polynomial_mod: ");
	print_polynomial(*negate_polynomial_mod(p, *differentiate(p)));
	printf("total_root_cnt: %d\n", total_root_cnt(p));
	printf("get_all_real_roots:\n");
	print_root_list(*get_all_real_roots(p, ROOT_2_ERR));
	printf("get_all_roots:\n");
	polynomial q = *calloc_polynomial(5);
	q.coefficients[0] = 1;
	q.coefficients[4] = -1;
	q.coefficients[5] = 1;
	print_complex_root_list(*get_all_roots(q, DEG5_ERR));
}

int main(int argc, char **argv) {
	test_root_functions();
	exit(0);
}

#endif
