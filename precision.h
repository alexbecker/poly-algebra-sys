// precision.h
// includes type declaration of entry_type and definitions of various precision constants

#ifndef PRECISION_H
#define PRECISION_H

// precision constants:
#define EVAL_ERR 1e-20
#define MIN_ROOT_ERR 1e-12
#define MATRIX_IND_ERR 1e-10

// the precision of some functions can be changed by changing the following typedefs:
typedef __float128 matrix_entry;
typedef double root_type;
typedef _Complex double complex;	// precision should be same as root_type

#endif
