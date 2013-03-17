// roots.h

#ifndef ROOTS_H
#define ROOTS_H

#include "polynomial.h"
#include "precision.h"

// accounts for error by treating root_type as a small ball
typedef struct ball {
	root_type center, radius;
} ball;

typedef struct complex_ball {
	complex center;
	root_type radius;
} complex_ball;

// array of roots of a polynomial
typedef struct root_list {
	int num_roots;
	ball *roots;	// if num_roots = 0, roots should be NULL
} root_list;

typedef struct complex_root_list {
	int num_roots;
	complex_ball *roots;
} complex_root_list;

ball *new_ball(root_type, root_type);

complex_ball *new_complex_ball(complex, root_type);

ball *copy_ball(ball);

root_list *root_to_root_list(root_type, root_type, int);

root_list *empty_root_list();

root_list *alloc_root_list(int);

complex_root_list *alloc_complex_root_list(int);

void free_root_list(root_list*);

void free_complex_root_list(complex_root_list*);

int test_ball_membership(root_type, ball);

int test_complex_ball_membership(complex, complex_ball);

root_list *merge_root_lists(root_list, root_list);

int test_for_root(polynomial, ball);

int root_cnt(polynomial, root_type, root_type);

int total_root_cnt(polynomial);

root_list *get_real_roots(polynomial, root_type, root_type, root_type);

root_list *get_all_real_roots(polynomial, root_type);

complex_root_list *get_all_roots(polynomial, root_type);

void print_ball(ball);

void print_complex_ball(complex_ball);

void print_root_list(root_list);

void print_complex_root_list(complex_root_list);

#endif
