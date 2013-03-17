// calculator.c
// stack-based calculator for manipulating algebraic numbers
// version: 0.2

#include "calc_interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
	printf("Welcome to algebraic calculator v. 0.2!\n"
			"Type \"help\" for a list of commands.\n");
	// create an empty stack
	stack s = (stack) malloc(sizeof(stack_entry));
	s->next = NULL;
	// perform user's operations indefinitely
	char command[25];
	while (1) {
		printf("Enter a command: ");
		// read input, skipping any newlines
		do {
			fgets(command, 25, stdin);
		} while (command[0] == '\n');
		if (!strcmp(command, "print\n")) {
			print(s);
		} else if (!strcmp(command, "print-full\n")) {
			print_full(s);
		} else if (!strcmp(command, "new\n")) {
			new(s);
		} else if (!strcmp(command, "refine\n")) {
			refine(s);
		} else if (!strcmp(command, "add\n")) {
			add(s);
		} else if (!strcmp(command, "negate\n")) {
			negate(s);
		} else if (!strcmp(command, "subtract\n")) {
			subtract(s);
		} else if (!strcmp(command, "multiply\n")) {
			multiply(s);
		} else if (!strcmp(command, "invert\n")) {
			invert(s);
		} else if (!strcmp(command, "divide\n")) {
			divide(s);
		} else if (!strcmp(command, "clear\n")) {
			clear(s);
		} else if (!strcmp(command, "clear-all\n")) {
			clear_all(s);
		} else if (!strcmp(command, "quit\n")) {
			exit(0);
		} else if (!strcmp(command, "help\n")) {
			help();
		} else {
			printf("Command not recognized. Type \"help\" for a list of commands.");
		}
	}
}
