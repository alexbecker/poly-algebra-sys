// algebraics.c
// implements basic functions for handling algebraic numbers

#include "algebraics.h"
#include "minpoly.h"
#include "resultant.h"
#include "factoring.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* ---------- Constructors ---------- */

static algebraic *ball_to_algebraic(ball b) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = *new_ball(b.center, b.radius);
	result->minimal_polynomial = NULL;
	
	return result;
}

// NOTE: this function requires user interaction
static algebraic *polynomial_to_algebraic(polynomial p, root_type error_bnd) {
	root_list *roots = get_all_real_roots(p, error_bnd);
	// only works if the polynomial has a root
	assert(roots->num_roots != 0);
	// ask the user which root they want to consider
	print_root_list(*roots);
	printf("Select a root [0-%d]: ", roots->num_roots - 1);
	int root_selected = roots->num_roots;
	while (root_selected >= roots->num_roots) {
		scanf("%d", &root_selected);
	}
	// create the algebraic
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = roots->roots[root_selected];
	result->minimal_polynomial = copy_polynomial(p);
	// free intermediates
	free_root_list(roots);
	
	return result;
}

// used if both the root and algebraic are known
static algebraic *both_to_algebraic(polynomial p, ball b) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = *new_ball(b.center, b.radius);
	result->minimal_polynomial = copy_polynomial(p);
	
	return result;
}

/* ---------- Memory Functions ---------- */

void free_algebraic(algebraic *a) {
	if (a->minimal_polynomial != NULL) {
		free_polynomial(a->minimal_polynomial);
	}
	free(a);
}

algebraic *copy_algebraic(algebraic a) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = *copy_ball(a.approx_val);
	if (a.minimal_polynomial == NULL) {
		result->minimal_polynomial = NULL;
	} else {
		result->minimal_polynomial = copy_polynomial(*(a.minimal_polynomial));
	}
	
	return result;
}

/* ---------- Testing / Modifying Functions ---------- */

// if the minimal polynomial is not already defined, compute it
// if get_min_poly does not find a minimal polynomial, it remains undefined
void define_minimal_polynomial(algebraic *a, int max_k, int max_deg) {
	if (a->minimal_polynomial == NULL)
		a->minimal_polynomial = get_min_poly(a->approx_val, max_k, max_deg);
}

// check whether an algebraic number is uniquely defined
// return 1 if so, 0 if not
// if the minimal polynomial is not defined, returns 1
int is_uniquely_defined(algebraic a) {
	if (a.minimal_polynomial != NULL) {
		// make sure there is exactly one root in the ball
		if (root_cnt(*(a.minimal_polynomial), a.approx_val.center - a.approx_val.radius, a.approx_val.center + a.approx_val.radius) != 1)
			return 0;
	}

	return 1;
}
		
// refine the approx_val for the minimal polynomial
// if multiple choices can be made, takes the smallest
// if the minimal polynomial is not defined, simply shrinks the radius of approx_val
void refine_approx_val(algebraic *a, root_type error) {
	if (a->approx_val.radius > error) {
		if (a->minimal_polynomial == NULL) {
			a->approx_val.radius = error;
		} else {
			root_list *refined_roots = get_real_roots(*(a->minimal_polynomial), a->approx_val.center - a->approx_val.radius, a->approx_val.center + a->approx_val.radius, error);
			// make sure a root exists
			assert(refined_roots->num_roots > 0);
			a->approx_val = refined_roots->roots[0];
			free_root_list(refined_roots);
		}
	}
}

/* ---------- Algebraic Arithmetic ---------- */

algebraic *add_algebraics(algebraic a, algebraic b) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = *new_ball(a.approx_val.center + b.approx_val.center, a.approx_val.radius + b.approx_val.radius);
	// if the minimal polynomials are not null, take their resultant
	if ((a.minimal_polynomial != NULL) && (b.minimal_polynomial != NULL)) {
		polynomial *min_poly_unfactored = resultant_sum(*(a.minimal_polynomial), *(b.minimal_polynomial));
		result->minimal_polynomial = find_factor(*min_poly_unfactored, result->approx_val);
		free_polynomial(min_poly_unfactored);
	} else {
		result->minimal_polynomial = NULL;
	}
	
	return result;
}

algebraic *negate_algebraic(algebraic a) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	result->approx_val = *new_ball(-a.approx_val.center, a.approx_val.radius);
	if (a.minimal_polynomial == NULL) {
		result->minimal_polynomial = NULL;
	} else {
		result->minimal_polynomial = linear_change_of_variables(*(a.minimal_polynomial), -1, 0);
	}
	
	return result;
}

algebraic *subtract_algebraics(algebraic a, algebraic b) {
	algebraic *b_negated = negate_algebraic(b);
	algebraic *result = add_algebraics(a, *b_negated);
	free_algebraic(b_negated);
	
	return result;
}

algebraic *mult_algebraics(algebraic a, algebraic b) {
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	// compute the new error
	root_type new_error = fabs(a.approx_val.center) * b.approx_val.radius + fabs(b.approx_val.center) * a.approx_val.radius;
	result->approx_val = *new_ball(a.approx_val.center * b.approx_val.center, new_error);
	// if the minimal polynomials are not null, take their resultant
	if ((a.minimal_polynomial != NULL) && (b.minimal_polynomial != NULL)) {
		polynomial *min_poly_unfactored = resultant_product(*(a.minimal_polynomial), *(b.minimal_polynomial));
		result->minimal_polynomial = find_factor(*min_poly_unfactored, result->approx_val);
		free_polynomial(min_poly_unfactored);
	} else {
		result->minimal_polynomial = NULL;
	}
	
	return result;
}

algebraic *invert_algebraic(algebraic a) {
	// test to make sure a is invertible
	assert(!test_ball_membership(0, a.approx_val));
	algebraic *result = (algebraic*) malloc(sizeof(algebraic));
	// compute the new error
	root_type new_error = a.approx_val.radius / (fabs(a.approx_val.center) - a.approx_val.radius);
	result->approx_val = *new_ball(1.0 / a.approx_val.center, new_error);
	if (a.minimal_polynomial == NULL) {
		result->minimal_polynomial = NULL;
	} else {
		result->minimal_polynomial = reverse_polynomial(*(a.minimal_polynomial));
	}
	
	return result;
}

algebraic *divide_algebraics(algebraic a, algebraic b) {
	algebraic *b_inverted = invert_algebraic(b);
	algebraic *result = mult_algebraics(a, *b_inverted);
	free_algebraic(b_inverted);
	
	return result;
}

/* ---------- Number-Theoretic Functions ---------- */

static complex_root_list *get_galois_conjugates(algebraic a) {
	assert(a.minimal_polynomial != NULL);
	return get_all_roots(*(a.minimal_polynomial), a.approx_val.radius);
}

/* ---------- Input / Output ---------- */

// reads a polynomial from stdin
// NOTE: this function requires user interaction
algebraic *read_algebraic() {
	algebraic *result;
	polynomial *entered_polynomial = NULL;
	double root, error;
	printf("Enter an approximate value for the algebraic number, or press enter if it is unknown: ");
	if (scanf("%lf", &root)) {	// the user has entered a number
		printf("Enter an upper bound for the error: ");
		scanf("%lf", &error);
		printf("Enter the minimal polynomial, or \"0\" if it is unknown: ");
		ball *entered_ball = new_ball(root, error);
		entered_polynomial = read_polynomial();
		if (entered_polynomial != NULL) {	// the user has entered a polynomial
			result = both_to_algebraic(*entered_polynomial, *entered_ball);
			free(entered_polynomial);
		} else {							// the user has not
			result = ball_to_algebraic(*entered_ball);
		}
		free(entered_ball);
	} else {					// the user has not
		// force the user to enter a polynomial
		while (entered_polynomial == NULL) {
			printf("Enter the minimal polynomial: ");
			entered_polynomial = read_polynomial();
		}
		printf("Enter an upper bound for the error of the root you want: ");
		scanf("%lf", &error);
		result = polynomial_to_algebraic(*entered_polynomial, error);
		free(entered_polynomial);
	}
	
	return result;
}

void print_algebraic(algebraic a) {
	print_ball(a.approx_val);
	if (a.minimal_polynomial == NULL) {
		printf("Minimal polynomial unknown.\n");
	} else {
		printf("Minimal polynomial: ");
		print_polynomial(*(a.minimal_polynomial));
	}
}

static void print_galois_conjugates(algebraic a) {
	printf("Galois conjugates:\n");
	complex_root_list *galois_conjugates = get_galois_conjugates(a);
	print_complex_root_list(*galois_conjugates);
	free_complex_root_list(galois_conjugates);
}

static void print_algebraic_degree(algebraic a) {
	if (a.minimal_polynomial == NULL) {
		printf("Degree unknown.\n");
	} else {
		printf("Degree: %d\n", a.minimal_polynomial->deg);
	}
}

void print_all_information(algebraic a) {
	print_algebraic(a);
	print_galois_conjugates(a);
	print_algebraic_degree(a);
}

/* ---------- Testing ---------- */
// to test, run "make test algebraics"

#ifdef TEST_ALGEBRAICS
#define DEG5_ALG -1.1673039782614186843
#define DEG5_ERR 1e-10

/* should output:
 * ball_to_algebraic:
 * Approximate value: -1.167304, Error: 0.000000
 * Minimal polynomial unknown.
 * define_minimal_polynomial:
 * Approximate value: -1.167304, Error: 0.000000
 * Minimal polynomial: x^5 - x^1 + 1
 * polynomial_to_algebraic:
 * Root 0: Approximate value: -1.414214, Error: 0.000000
 * Root 1: Approximate value: 1.414214, Error: 0.000000
 * Approximate value: -1.414214, Error: 0.000000
 * Minimal polynomial: x^2 - 2
 * add_algebraics:
 * Approximate value: -2.581518, Error: 0.000000
 * Minimal polynomial: x^10 - 10x^8 + 38x^6 + 2x^5 - 100x^4 + 40x^3 + 121x^2 + 38x^1 - 17
 * subtract_algebraics:
 * Approximate value: 0.246910, Error: 0.000000
 * Minimal polynomial: x^10 - 10x^8 + 38x^6 + 2x^5 - 100x^4 + 40x^3 + 121x^2 + 38x^1 - 17
 * mult_algebraics:
 * Approximate value: 1.650817, Error: 0.000000
 * Minimal polynomial: 2107536x^10 - 1022x^9 + 24532x^8 - 337314x^7 - 13927955x^6 - 16764517x^5 + 63502981x^4 - 156288926x^3 + 270743008x^2 - 197531212x^1
 * divide_algebraics:
 * Approximate value: 0.825409, Error: 0.000000
 * Minimal polynomial: -1251375499x^10 + 593x^9 - 14225x^8 + 195592x^7 + 623987515x^6 + 9720926x^5 - 36822283x^4 + 90624330x^3 - 215648899x^2 + 114538721x^1
 * print_galois_conjugates:
 * Galois conjugates:
 * Root 0: Approximate value: -1.414214, Error: 0.000000
 * Root 1: Approximate value: 1.414214, Error: 0.000000
 * read_algebraic can be checked by hand */
void test_algebraics_functions() {
	printf("ball_to_algebraic:\n");
	algebraic *a = ball_to_algebraic(*new_ball(DEG5_ALG, DEG5_ERR));
	print_algebraic(*a);
	printf("define_minimal_polynomial:\n");
	define_minimal_polynomial(a, 3, 5);
	print_algebraic(*a);
	polynomial *p = calloc_polynomial(2);
	p->coefficients[0] = 1;
	p->coefficients[2] = -2;
	// NOTE: "should output" expects you to enter 0
	printf("polynomial_to_algebraic:\n");
	algebraic *b = polynomial_to_algebraic(*p, 1e-10);
	print_algebraic(*b);
	printf("add_algebraics:\n");
	print_algebraic(*add_algebraics(*a,*b));
	printf("subtract_algebraics:\n");
	print_algebraic(*subtract_algebraics(*a,*b));
	printf("mult_algebraics:\n");
	print_algebraic(*mult_algebraics(*a,*b));
	printf("divide_algebraics:\n");
	print_algebraic(*divide_algebraics(*a,*b));
	printf("print_galois_conjugates:\n");
	print_galois_conjugates(*b);
	printf("read_algebraic:\n");
	print_algebraic(*read_algebraic());
}

int main(int argc, char **argv) {
	test_algebraics_functions();
	exit(0);
}

#endif
