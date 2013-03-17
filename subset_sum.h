#ifndef SUBSET_SUM_H
#define SUBSET_SUM_H

// get declaration of entry_type
#include "precision.h"

typedef struct subset_sum_problem {
	int size;
	root_type *items;
	root_type sum;
} subset_sum_problem;

int *subset_sum_certificate(subset_sum_problem, root_type *);

#endif
