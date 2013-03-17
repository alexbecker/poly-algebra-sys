// interface.h

#ifndef CALC_INTERFACE_H
#define CALC_INTERFACE_H

#include "algebraics.h"

typedef struct stack_entry {
	algebraic *value;
	struct stack_entry *next;	// NULL at bottom of stack
} stack_entry, *stack;	// stack should always be a pointer to the top

void print(stack);

void print_full(stack);

void new(stack);

void refine(stack);

void add(stack);

void negate(stack);

void subtract(stack);

void multiply(stack);

void invert(stack);

void divide(stack);

void clear(stack);

void clear_all(stack);

void help();

#endif
