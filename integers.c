// integers.c
// impliments basic integer functions
// must be compiled according to C99 for % to behave appropriately

#include "integers.h"
#include <stdlib.h>
#include <stdio.h>

int gcd(int a, int b) {
	if (a % b == 0) {
		return b;
	} else if (b % a == 0) {
		return a;
	} else if (a == 1 || a == -1 || b == 1 || b == -1) {
		return 1;
	} else {
		return gcd(b, a % b);	// orders switched in case b > a
	}
}

// used only for testing
// to test, run "make test integers"
#ifdef TEST_INTEGERS

// should output 3, 4 up to sign
void test_gcd() {
	printf("gcd: %d, %d\n", gcd(12, 27), gcd(80, -12));
}

int main(int argc, char **argv) {
	test_gcd();
	exit(0);
}

#endif
