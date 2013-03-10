// rationals.c
// defines rationals and impliments basic functions

#include "rationals.h"
#include "integers.h"
#include <stdlib.h>
#include <stdio.h>

rational *to_rational(int num, int denom) {
	rational *result = (rational*) malloc(sizeof(rational));
	result->num = num;
	resutl->denom = denom;
	
	return result;
}

// casts a rational to an integer, assuming it equals an integer
int to_int(rational r) {
	reduce(r);	// reduce to lowest terms
	assert(r.denom == 1);	// assert that the rational equals an integer
	
	return r.num;
}

// reduces a rational to lowest terms and makes the denominator positive
void reduce(rational *r) {
	int common_factor = gcd(r->num, r->denom);
	r->num /= common_factor;
	r->denom /= common_factor;
	if (r->denom < 0) {
		r->num *= -1;
		r->denom *= -1;
	}
}

// prints a rational to stdout in human-readable form
void print_rational(rational r) {
	printf("%d/%d", r.num, r.denom);
}

#ifdef TEST_RATIONALS

void test_reduce() {
	rational r = to_rational(8, -12);
	reduce(&r);
	print_rational(r);
}

int main(int argc, char ** argv) {
	test_reduce();
	exit(0);
}

#endif
