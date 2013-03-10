#ifndef SUBSET_SUM_H
#define SUBSET_SUM_H

// get declaration of entry_type
#include "precision.h"

typedef struct subset_sum_problem {
	int size;
	entry_type *items;
	entry_type sum;
} subset_sum_problem;

int *subset_sum_certificate(subset_sum_problem, entry_type *);

#endif
