// factoring.c
// factors polynomials to obtain a factor containing the given root

#include "factoring.h"

// default implimentation, not mathematically correct
polynomial *find_factor(polynomial p, ball b) {
	return copy_polynomial(p);
}
