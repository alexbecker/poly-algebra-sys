// subset_sum.c
// solves the subset sum problem in O(2^(n/2)) time and space

#include "subset_sum.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>

// merges two sorted lists of size n in O(n) time
// breaks ties towards first list
// frees its inputs
// NOTE: list1 and list2 must have memory for an additional entry
// the returned list has memory for an additional memory
static entry_type *merge_lists(int size, entry_type *list1, entry_type *list2) {
	entry_type *result = (entry_type*) malloc(sizeof(entry_type) * (2 * size + 1));
	assert(result != NULL);
	
	// add huge last entries to list1 and list2 to avoid checking for the end of the list
	list1[size] = (entry_type) DBL_MAX / 2.0;
	list1[size] = (entry_type) DBL_MAX / 2.0;
	// keep track of position in list 1 and 2
	int list1_ind = 0, list2_ind = 0;
	for (int i=0; i<2*size+1; i++) {
		if (list1[list1_ind] <= list2[list2_ind]) {
			result[i] = list1[list1_ind++];
		} else {
			result[i] = list2[list2_ind++];
		}
	}
	
	free(list1);
	free(list2);
	return result;
}

// used to get list of sorted sums of the two subsets
static entry_type *get_sorted_sums(int size, entry_type *vals) {
	int num_sums = 1;
	entry_type *list1 = (entry_type*) malloc(sizeof(entry_type) * (num_sums + 1));
	list1[0] = 0;
	for (int i=0; i<size; i++) {
		entry_type *list2 = (entry_type*) malloc(sizeof(entry_type) * (num_sums + 1));
		for (int j=0; j<num_sums; j++)
			list2[j] = list1[j] + vals[i];
		list1 = merge_lists(num_sums, list1, list2);
		num_sums *= 2;
	}
	return list1;
}

// returns smallest absolute difference between given sum and a subset sum
// runs in O(2^(n/2)) time and space
static double subset_sum(int size, entry_type *vals, entry_type sum) {
	// the set is split into two equal-size subsets
	int size1 = size / 2;
	int size2 = size - size1;
	entry_type *vals1 = (entry_type*) malloc(sizeof(entry_type) * size1);
	entry_type *vals2 = (entry_type*) malloc(sizeof(entry_type) * size2);
	for (int i=0; i<size1; i++)
		vals1[i] = vals[i];
	for (int i=0; i<size2; i++)
		vals2[i] = vals[i + size1];
	
	entry_type *sorted_sums1 = get_sorted_sums(size1, vals1);
	entry_type *sorted_sums2 = get_sorted_sums(size2, vals2);
	int num_sorted_sums1 = (int) round(pow(2.0, (double) size1));
	int num_sorted_sums2 = (int) round(pow(2.0, (double) size2));
	// vals no longer needed, so are freed
	free(vals1);
	free(vals2);
	
	int ind1 = 0, ind2 = num_sorted_sums2 - 1;
	double smallest_err = DBL_MAX;
	
	while (1) {
		entry_type cur_dif = sum - (sorted_sums1[ind1] + sorted_sums2[ind2]);
		if (fabs((double) cur_dif) < smallest_err)
			smallest_err = fabs((double) cur_dif);
		if (cur_dif > (entry_type) 0) {	// subset sum is too low
			if (ind1 < num_sorted_sums1 - 1) {
				ind1++;	// move to larger element of sorted_sums1
			} else {
				break;	// finished searching
			}
		} else if (cur_dif < (entry_type) 0) {	// subset sum is too high
			if (ind2 > 0) {
				ind2--;	// move to smaller element of sorted_sums2
			} else {
				break;	// finished searching
			}
		} else {	// subset sum is exactly given sum
			break;
		}
	}
	
	free(sorted_sums1);
	free(sorted_sums2);
	return smallest_err;
}

// returns the subset which sums closest to error, in the form of a list of 0s and 1s indicating inclusion and exclusion
// records the error in the passed error parameter
// runs in O(2^((n+1)/2)) time and O(2^(n/2)) space
int *subset_sum_certificate(subset_sum_problem problem, entry_type *error) {
	// initialize the set of indices being checked to all indices
	int *inds = (int*) malloc(sizeof(int) * problem.size);
	for (int i=0; i<problem.size; i++)
		inds[i] = 1;
	
	// the minimum acheivable error is that found searching the whole set
	*error = subset_sum(problem.size, problem.items, problem.sum);
	
	// search through the tree of subsets, selecting the appropriate branch at each step using subset_sum
	entry_type *subset = malloc(sizeof(entry_type) * problem.size);
	int subset_size;
	for (int i=0; i<problem.size; i++) {
		subset_size = 0;
		// first try leaving out the ith element
		inds[i] = 0;
		for (int j=0; j<problem.size; j++) {
			if (inds[j] == 1)
				subset[subset_size++] = problem.items[j];
		}
		// test if we can still get within error of sum
		if (subset_sum(subset_size, subset, problem.sum) > *error + EVAL_ERR)
			inds[i] = 1;	// if not, puts the ith element back in
	}
	
	free(subset);
	return inds;
}

// used only for testing
// to test, run "make test subset_sum"
#ifdef TEST_SUBSET_SUM

// should output 0-19 in order
void test_merge_lists() {
	entry_type *list1 = (entry_type*) malloc(sizeof(entry_type) * 11);
	entry_type *list2 = (entry_type*) malloc(sizeof(entry_type) * 11);
	for (int i=0; i<10; i++) {
		list1[i] = 2*i;
		list2[i] = 2*i+1;
	}
	entry_type *list3 = merge_lists(10, list1, list2);
	printf("merge_lists: ");
	for (int i=0; i<20; i++) {
		printf("%lf\n", (double) list3[i]);
	}
}

// should output 0-15 in order
void test_get_sorted_sums() {
	entry_type *list = (entry_type*) malloc(sizeof(entry_type) * 5);
	list[0] = 1;
	for (int i=1; i<4; i++) {
		list[i] = 2*list[i-1];
	}
	entry_type *sorted_sums = get_sorted_sums(4, list);
	printf("subset_sum: ");
	for (int i=0; i<16; i++) {
		printf("%lf\n", (double) sorted_sums[i]);
	}
}

// should output 0.1
void test_subset_sum() {
	entry_type *list = (entry_type*) malloc(sizeof(entry_type) * 10);
	list[0] = 1;
	for (int i=1; i<10; i++) {
		list[i] = 2*list[i-1];
	}
	printf("subset_sum: %lf\n", subset_sum(10, list, (entry_type) 10.1));
}

// should output error of .1 and subset 0,0,1,0,0,1,1,0,0,0 
void test_subset_sum_certificate() {
	entry_type *list = (entry_type*) malloc(sizeof(entry_type) * 10);
	list[0] = 1;
	for (int i=1; i<10; i++) {
		list[i] = 2*list[i-1];
	}
	entry_type error;
	subset_sum_problem *problem = (subset_sum_problem*) malloc(sizeof(subset_sum_problem));
	problem->size = 10;
	problem->items = list;
	problem->sum = (entry_type) 100.1;
	int *subset = subset_sum_certificate(*problem, &error);
	printf("subset_sum_certificate: Error: %lf\nSubset: %d", (double) error, subset[0]);
	for (int i=1; i<10; i++)
		printf(", %d", subset[i]);
	printf("\n");
}

// used to test all functions
int main(int argc, char **argv) {
	test_merge_lists();
	test_get_sorted_sums();
	test_subset_sum();
	test_subset_sum_certificate();
	exit(0);
}

#endif
