# Makefile for minimal polynomial project

CC = gcc
CFLAGS = -std=gnu99
LDFLAGS = -lm
SRC = minpoly.c subset_sum.c polynomial.c matrices.c integers.c rationals.c
OBJ = minpoly.o subset_sum.o polynomial.o matrices.o integers.o rationals.o

minpoly: CFLAGS += -Wall -DTEST_MINPOLY
minpoly: minpoly.o subset_sum.o polynomial.o

subset_sum: CFLAGS += -Wall -DTEST_SUBSET_SUM
subset_sum: subset_sum.o

polynomial: CFLAGS += -Wall -DTEST_POLYNOMIAL
polynomial: polynomial.o

matrices: CFLAGS += -Wall -DTEST_MATRICES
matrices: matrices.o

integers: CFLAGS += -Wall -DTEST_INTEGERS
integers: integers.o

rationals: CFLAGS += Wall -DTEST_RATIONALS
rationals: rationals.o

# remove object files prior to compiling test versions
test:
	rm -f ${OBJ}

depend:
	./remove_depends.sh
	gcc -MM ${SRC} >> Makefile

clean:
	rm -f ${OBJ} minpoly subset_sum polynomial

# dependencies listed by gcc -MM
minpoly.o: minpoly.c minpoly.h polynomial.h precision.h subset_sum.h
subset_sum.o: subset_sum.c subset_sum.h precision.h
polynomial.o: polynomial.c polynomial.h precision.h
matrices.o: matrices.c matrices.h precision.h
integers.o: integers.c integers.h
