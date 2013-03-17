// calc_interface.c
// interface between calculator.c and algebraics.c

#include "calc_interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* ---------- Stack Functions ---------- */

static void push(algebraic *a, stack s) {
	stack s_cpy = (stack) malloc(sizeof(stack_entry));
	s_cpy->value = s->value;
	s_cpy->next = s->next;
	s->value = a;
	s->next = s_cpy;
}

static algebraic *pop(stack s) {
	assert(s->value != NULL);	// assert that stack is not empty
	algebraic *result = s->value;
	if (s->next != NULL) {
		s->value = s->next->value;
		free(s->next);
		s->next = s->next->next;	// hanging pointer
	} else {
		s->value = NULL;
	}

	return result;
}

/* ---------- Operations ---------- */

void print(stack s) {
	algebraic *top = pop(s);
	print_algebraic(*top);
	push(top, s);
}

void print_full(stack s) {
	algebraic *top = pop(s);
	print_all_information(*top);
	push(top, s);
}

void new(stack s) {
	algebraic *new_algebraic = read_algebraic();
	push(new_algebraic, s);
}

void refine(stack s) {
	printf("Would you like to improve the accuracy of the algebraic number, or find its minimal polynomial? [0-1]: ");
	int selection = 2;
	while (selection > 1) {
		scanf("%d", &selection);
	}
	if (selection == 0) {
		double new_error;
		printf("Enter the new error: ");
		scanf("%lf", &new_error);
		refine_approx_val(s->value, new_error);
	} else {
		int max_k, max_n;
		printf("Enter the maximum exponent on the coefficients to consider: ");
		scanf("%d", &max_k);
		printf("Enter the maximum degree to consider: ");
		scanf("%d", &max_n);
		define_minimal_polynomial(s->value, max_k, max_n);
	}
}

void add(stack s) {
	algebraic *a = pop(s), *b = pop(s);
	push(add_algebraics(*a, *b), s);
	free_algebraic(a);
	free_algebraic(b);
}

void negate(stack s) {
	algebraic *a = pop(s);
	push(negate_algebraic(*a), s);
	free_algebraic(a);
}

void subtract(stack s) {
	algebraic *a = pop(s), *b = pop(s);
	push(subtract_algebraics(*a, *b), s);
	free_algebraic(a);
	free_algebraic(b);
}

void multiply(stack s) {
	algebraic *a = pop(s), *b = pop(s);
	push(mult_algebraics(*a, *b), s);
	free_algebraic(a);
	free_algebraic(b);
}

void invert(stack s) {
	algebraic *a = pop(s);
	push(invert_algebraic(*a), s);
	free_algebraic(a);
}

void divide(stack s) {
	algebraic *a = pop(s), *b = pop(s);
	push(divide_algebraics(*a, *b), s);
	free_algebraic(a);
	free_algebraic(b);
}

void clear(stack s) {
	if (s->next != NULL) {
		algebraic *a = pop(s);
		free_algebraic(a);
	}
}

void clear_all(stack s) {
	while (s->next != NULL) {
		algebraic *a = pop(s);
		free_algebraic(a);
	}
}

void help() {
	printf("The following commands are available:\n"
			"print - prints the value on top of the stack\n"
			"print-full - prints all known information about the value on top of the stack\n"
			"new - allows the user to define a new algebraic number on top of the stack\n"
			"refine - refine the value on top of the stack in one of several ways\n"
			"add - adds the two values on top of the stack\n"
			"negate - negates the value on top of the stack\n"
			"subtract - subtracts the two values on top of the stack\n"
			"multiply - multiplies the two values on top of the stack\n"
			"invert - inverts the value on top of the stack\n"
			"divide - divides the two values on top of the stack\n"
			"clear - deletes the value on top of the stack\n"
			"clear-all - deletes all values on the stack\n"
			"quit - exits the program\n"
			"help - prints this message\n");
}
