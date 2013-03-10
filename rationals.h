#ifndef RATIONALS_H
#define RATIONALS_h

typedef struct rational {
	int num, denom;
} rational;

rational *to_rational(int);

int to_int(rational);

void reduce(rational*);

void print_rational(rational);

#endif 
