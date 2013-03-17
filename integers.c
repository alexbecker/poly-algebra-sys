// integers.c
// implements basic integer functions
// must be compiled according to C99 for % to behave appropriately

#include "integers.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#define CONT_FRAC_LEN 15

int gcd(int a, int b) {
	if (b == 0)
		return a;
	return gcd(b, a % b);	// orders switched in case b > a
}

int lcm(int a, int b) {
	if (a == 0)
		return 0;
	return (a * b) / gcd(a, b);
}

// helper function for denom
static matrix_entry eval_continued_fraction(int *fraction, int length) {
	if (length == 1)
		return (matrix_entry) fraction[0];
	return fraction[0] + 1 / eval_continued_fraction(fraction + 1, length - 1);
}

// returns the denominator of d as a fraction
// computes successive continued fraction approximations until one is close enough
int denom(matrix_entry d) {
/*	if (d < 0)	// handle negative d
		return denom(-d);
*/	// initialize the continued fraction
	int length = 0, max_length = CONT_FRAC_LEN;
	int *fraction = (int*) malloc(sizeof(int) * max_length);
	fraction[0] = (int) floor(d);
	matrix_entry next_term_unrounded = 1 / (d - fraction[0]);
	// increase the precesion of the continued fraction
	while (fabsl(d - eval_continued_fraction(fraction, ++length)) > MATRIX_IND_ERR) {
		fraction[length] = (int) next_term_unrounded;
		next_term_unrounded = 1.0 / (next_term_unrounded - fraction[length]);
		// if necessary, allocate more memory for the continued fraction
		if (length == max_length) {
			max_length *= 2;
			fraction = realloc(fraction, sizeof(int) * max_length);
		}
	}
	// get the denominator of the continued fraction
	int result;
	if (length == 1) {
		result = 1;
	} else {	// uses the fundamental recurrence formulas
		int *denoms = (int*) malloc(sizeof(int) * 2);
		denoms[0] = 1;
		denoms[1] = fraction[1];
		for (int i=2; i<length; i++) {
			denoms[i % 2] = fraction[i] * denoms[(i - 1) % 2] + denoms[(i - 2) % 2];
		}
		result = denoms[(length - 1) % 2];
		free(denoms);
	}
	free(fraction);
	
	return result;
}

// to test, run "make test integers"
#ifdef TEST_INTEGERS

// should output 3, -4
void test_gcd() {
	printf("gcd: %d, %d\n", gcd(12, 27), gcd(80, -12));
}

// should output 220, -36
void test_lcm() {
	printf("lcm: %d, %d\n", lcm(55, 4), lcm(-9, 12));
}

// should output 187, 781
void test_denom() {
	printf("denom: %d, %d\n", denom(3.0 / 561.0), denom(-5571.0 / 781.0));
}

int main(int argc, char **argv) {
	test_gcd();
	test_lcm();
	test_denom();
	exit(0);
}

#endif
