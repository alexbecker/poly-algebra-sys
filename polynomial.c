// polynomial.c
// implements basic functions for handling integer polynomials

#include "polynomial.h"
#include "integers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#define abs(a) ((a >= 0) ? a : -a)	// only use for integers

/* ---------- Constructors ---------- */

// casts an int to a polynomial of degree 0
polynomial *int_to_polynomial(int n) {
	polynomial *result = alloc_polynomial(0);
	result->coefficients[0] = n;
	
	return result;
}

/* ---------- Memory Functions ---------- */

// allocates a polynomial on the heap
polynomial *alloc_polynomial(int degree) {
	polynomial *result = (polynomial*) malloc(sizeof(polynomial));
	result->deg = degree;
	result->coefficients = (int*) malloc(sizeof(int) * (degree + 1));
	
	return result;
}

// allocates a polynomial on the heap with all coefficients 0
polynomial *calloc_polynomial(int degree) {
	polynomial *result = alloc_polynomial(degree);
	for (int i=0; i<=degree; i++) {
		result->coefficients[i] = 0;
	}
	
	return result;
}

// given a polynomial p, makes a copy of p on the heap
polynomial *copy_polynomial(polynomial p) {
	polynomial *p_cpy = alloc_polynomial(p.deg);
	for (int i=0; i<=p.deg; i++) {
		p_cpy->coefficients[i] = p.coefficients[i];
	}
	
	return p_cpy;
}

// frees a polynomial
void free_polynomial(polynomial *p) {
	free(p->coefficients);
	free(p);
}

// strips leading zeros from a polynomial
// many functions require that a polynomial have no leading zeros
void strip_leading_zeros(polynomial *p) {
	// determine the number of leading zeros
	int num_leading_zeros = 0;
	while (p->coefficients[num_leading_zeros] == 0) {
		num_leading_zeros++;
	}
	// check is the polynomial was all zeros
	if (num_leading_zeros > p->deg) {
		p->deg = 0;
	} else {
		int *new_coefficients = (int*) malloc(sizeof(int) * (p->deg - num_leading_zeros + 1));
		for (int i=0; i<=p->deg-num_leading_zeros; i++) {
			new_coefficients[i] = p->coefficients[i + num_leading_zeros];
		}
		p->deg -= num_leading_zeros;
		free(p->coefficients);
		p->coefficients = new_coefficients;
	}
}

/* ---------- Polynomial Arithmetic / Calculus ---------- */

// evaluate a polynomial using horner's scheme
root_type eval_polynomial(polynomial p, root_type x) {
	if (p.deg == 0)
		return p.coefficients[0];
	p.deg--;	// effectively decrease degree of p by 1
	
	return p.coefficients[p.deg + 1] + x * eval_polynomial(p, x);
}

// complex polynomial evaluation
complex eval_polynomial_complex(polynomial p, complex x) {
	if (p.deg == 0)
		return p.coefficients[0];
	p.deg--;	// effectively decrease degree of p by 1
	
	return p.coefficients[p.deg + 1] + x * eval_polynomial_complex(p, x);
}

// adds two polynomials
polynomial *add_polynomials(polynomial p, polynomial q) {
	polynomial *result;
	if (p.deg > q.deg) {
		result = copy_polynomial(p);
		for (int i=p.deg-q.deg; i<=p.deg; i++) {
			result->coefficients[i] += q.coefficients[i - p.deg + q.deg];
		}
	} else {
		result = copy_polynomial(q);
		for (int i=q.deg-p.deg; i<=q.deg; i++) {
			result->coefficients[i] += p.coefficients[i - q.deg + p.deg];
		}
	}
	strip_leading_zeros(result);
	
	return result;
}

polynomial *negate_polynomial(polynomial p) {
	polynomial *result = copy_polynomial(p);
	for (int i=0; i<=p.deg; i++) {
		result->coefficients[i] *= -1;
	}
	
	return result;
}

polynomial *subtract_polynomials(polynomial p, polynomial q) {
	polynomial *q_neg = negate_polynomial(q);
	polynomial *result = add_polynomials(p, *q_neg);
	free_polynomial(q_neg);
	
	return result;
}

// returns p with each coefficient mulitplied by n
polynomial *scalar_mult_polynomial(polynomial p, int n) {
	polynomial *result = copy_polynomial(p);
	for (int i=0; i<=p.deg; i++) {
		result->coefficients[i] *= n;
	}
	
	return result;
}

// returns p*x^n
polynomial *increase_degree(polynomial p, int n) {
	polynomial *result = calloc_polynomial(p.deg + n);
	for(int i=0; i<=p.deg; i++) {
		result->coefficients[i] = p.coefficients[i];
	}
	
	return result;
}

// returns x^p.deg*p(1/x)
polynomial *reverse_polynomial(polynomial p) {
	polynomial *result = alloc_polynomial(p.deg);
	for (int i=0; i<=p.deg; i++) {
		result->coefficients[i] = p.coefficients[p.deg - i];
	}
	strip_leading_zeros(result);
	
	return result;
}

polynomial *mult_polynomials(polynomial p, polynomial q) {
	polynomial *result = calloc_polynomial(p.deg + q.deg);
	for (int i=0; i<=result->deg; i++) {
		for (int j=0; j<=i; j++) {
			if ((j <= p.deg) && (i - j <= q.deg))
				result->coefficients[result->deg - i] += p.coefficients[p.deg - j] * q.coefficients[q.deg - (i - j)];
		}
	}
	
	return result;
}

// compute powers of p by repeated squaring
polynomial *polynomial_power(polynomial p, int n) {
	// base cases
	if (n == 0)
		return int_to_polynomial(1);
	if (n == 1)
		return copy_polynomial(p);
	// inductive case
	polynomial *half_pow = polynomial_power(p, n / 2);
	polynomial *result, *temp_result;	// temp_result holds p^(n-1)
	if (n % 2 == 0) {
		result = mult_polynomials(*half_pow, *half_pow);
	} else {
		temp_result = mult_polynomials(*half_pow, *half_pow);
		result = mult_polynomials(*temp_result, p);
		free_polynomial(temp_result);
	}
	free_polynomial(half_pow);
	
	return result;
}

polynomial *compose_polynomials(polynomial p, polynomial q) {
	polynomial *result = calloc_polynomial(p.deg * q.deg), *q_power = int_to_polynomial(1);
	for (int i=0; i<=p.deg; i++) {
		polynomial *scaled_q_power = scalar_mult_polynomial(*q_power, p.coefficients[p.deg - i]);
		polynomial *temp_result = add_polynomials(*result, *scaled_q_power);
		free_polynomial(scaled_q_power);
		free_polynomial(result);
		result = temp_result;
		polynomial *temp_q_power = mult_polynomials(*q_power, q);
		free_polynomial(q_power);
		q_power = temp_q_power;
	}
	
	return result;
}

// given p, a, b returns p(ax+b)
polynomial *linear_change_of_variables(polynomial p, int a, int b) {
	polynomial *result;
	if (a == 0) {	// deal with case a = 0
		result = alloc_polynomial(0);
		result->coefficients[0] = b;
	} else {	// other cases
		polynomial *lin_poly = alloc_polynomial(1);
		lin_poly->coefficients[0] = a;
		lin_poly->coefficients[1] = b;
		result = compose_polynomials(p, *lin_poly);
		free_polynomial(lin_poly);
	}
	
	return result;
}
polynomial *differentiate(polynomial p) {
	if (p.deg == 0)	// take care of constant polynomials
		return int_to_polynomial(0);
	polynomial *p_prime = alloc_polynomial(p.deg - 1);
	for (int i=0; i<p.deg; i++) {
		p_prime->coefficients[i] = (p.deg - i) * p.coefficients[i];
	}
	
	return p_prime;
}

// given polynomials p and q, returns the remainder of p mod q
// the result is not the true remainder, but rather the true remainder with cleared denominators
polynomial *polynomial_mod(polynomial p, polynomial q) {
	// check if there is anything to do
	if ((p.deg < q.deg) || (p.deg == 0 && p.coefficients[0] == 0))
		return copy_polynomial(p);	// if not, result is identical to p
	
	int pq_gcd = gcd(p.coefficients[0], q.coefficients[0]);
	polynomial *p_modified = scalar_mult_polynomial(p, q.coefficients[0] / pq_gcd);
	polynomial *q_modified = scalar_mult_polynomial(q, p.coefficients[0] / pq_gcd);
	polynomial *q_increased_deg = increase_degree(*q_modified, p.deg - q.deg);
	polynomial *temp_result = subtract_polynomials(*p_modified, *q_increased_deg);
	free_polynomial(p_modified);
	free_polynomial(q_modified);
	free_polynomial(q_increased_deg);
	polynomial *result = polynomial_mod(*temp_result, q);
	free_polynomial(temp_result);
	
	return result;
}

/* ---------- Input / Output ---------- */

// reads a polynomial from stdin
// returns NULL on failure
// NOTE: requires terms be in order of decresing degree
// NOTE: requires highest order coefficient is nonzero
// NOTE: every nonzero term must be prefaced and followed by a number, including 1
// NOTE: terms cannot be seperated by -, e.g. use 2x^5 + -1x^1 instead of 2x^5 - 1x^1
// NOTE: 0 is considered a null polynomial
// SUPPORTED VARIABLES: xyzw
polynomial *read_polynomial() {
	char input[256], *term_or_deg;
	// read input, skipping any newlines
	do {
		fgets(input, 256, stdin);
	} while (input[0] == '\n');
	// check if a polynomial was entered, if not return NULL
	if (strlen(input) <= 1)
		return NULL;
	// check if 0 was entered, if so return NULL
	if (input[0] == '0')
		return NULL;
	// scan through the input string, using strtok to break up the string
	int max_deg, degree, coeff;
	// read the first coefficient, if unsucessful return NULL
	term_or_deg = strtok(input, "xyzw^ +");
	if (!sscanf(term_or_deg, "%d", &coeff))
		return NULL;
	// read the highest degree
	term_or_deg = strtok(NULL, "xyzw^ +");
	if (!sscanf(term_or_deg, "%d", &max_deg))
		return NULL;
	// allocate the polynomial based on this information
	polynomial *result = calloc_polynomial(max_deg);
	result->coefficients[0] = coeff;
	// scan through the rest of the input
	while ((term_or_deg = strtok(NULL, "xyzw^ +")) != NULL) {
		// read the coefficient
		if (!sscanf(term_or_deg, "%d", &coeff))
			return NULL;
		// read the degree
		term_or_deg = strtok(NULL, "xyzw^ +");
		// if NULL, the degree is 0
		if (term_or_deg == NULL) {
			degree = 0;
		} else {
			if (!sscanf(term_or_deg, "%d", &degree))
				return NULL;
		}
		// set the corresponding coefficient
		result->coefficients[max_deg - degree] = coeff;
	}
	
	return result;
}

// prints a polynomial to stdout in human-readable form
void print_polynomial(polynomial p) {
	if (p.deg == 0) {	// special case of degree 0 polynomial
		printf("%d\n", p.coefficients[0]);
		return;
	}
	if (p.coefficients[0] == 1) {
		printf("x^%d", p.deg);
	} else if (p.coefficients[0] == -1) {
		printf("-x^%d", p.deg);
	} else {
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

/* ---------- Testing ---------- */
// to test, run "make test polynomial"

#ifdef TEST_POLYNOMIAL
#define ROOT_2 1.414213562373095

/* should output:
 * print_polynomial: x^2 - 2
 * stript_leading_zeros: x^1 + 1
 * eval_polynomial: 0
 * add_polynomials: x^3 + x^2 - 3
 * negate_polynomial: -x^2 + 2
 * subtract_polynomials: -x^3 + x^2 - 1
 * scalar_mult_polynomials: 3x^2 - 6
 * increase_degree: x^4 - 2x^2
 * reverse_polynomial: -2x^2 + 1
 * mult_polynomials: x^5 - 2x^3 - x^2 + 2
 * polynomial_powers: x^6 - 6x^4 + 12x^2 - 8
 * compose_polynomials: x^6 - 2x^3 - 1
 * linear_change_of_variables: x^2 - 2x^1 - 1
 * differentiate: 2x^1
 * polynomial_mod: 2x^1 - 1	
 * read_polynomial can be checked by hand */
void test_polynomial_functions() {
	polynomial *p = alloc_polynomial(2);
	p->coefficients[0] = 1;
	p->coefficients[1] = 0;
	p->coefficients[2] = -2;
	polynomial *q = alloc_polynomial(3);
	q->coefficients[0] = 1;
	q->coefficients[1] = 0;
	q->coefficients[2] = 0;
	q->coefficients[3] = -1;
	printf("print_polynomial: ");
	print_polynomial(*p);
	polynomial *malformed = alloc_polynomial(2);
	malformed->coefficients[0] = 0;	// leading coefficient is 0
	malformed->coefficients[1] = 1;
	malformed->coefficients[2] = 1;
	strip_leading_zeros(malformed);
	printf("stript_leading_zeros: ");
	print_polynomial(*malformed);
	printf("eval_polynomial: %lf\n", (double) eval_polynomial(*p, ROOT_2));
	printf("add_polynomials: ");
	print_polynomial(*add_polynomials(*p, *q));
	printf("negate_polynomial: ");
	print_polynomial(*negate_polynomial(*p));
	printf("subtract_polynomials: ");
	print_polynomial(*subtract_polynomials(*p, *q));
	printf("scalar_mult_polynomials: ");
	print_polynomial(*scalar_mult_polynomial(*p, 3));
	printf("increase_degree: ");
	print_polynomial(*increase_degree(*p, 2));
	printf("reverse_polynomial: ");
	print_polynomial(*reverse_polynomial(*p));
	printf("mult_polynomials: ");
	print_polynomial(*mult_polynomials(*p, *q));
	printf("polynomial_powers: ");
	print_polynomial(*polynomial_power(*p, 3));
	printf("compose_polynomials: ");
	print_polynomial(*compose_polynomials(*p, *q));
	printf("linear_change_of_variables: ");
	print_polynomial(*linear_change_of_variables(*p, -1, 1));
	printf("differentiate: ");
	print_polynomial(*differentiate(*p));
	printf("polynomial_mod: ");
	print_polynomial(*polynomial_mod(*q, *p));
	printf("read_polynomial: ");
	print_polynomial(*read_polynomial());
}

int main(int argc, char **argv) {
	test_polynomial_functions();
	exit(0);
}

#endif
