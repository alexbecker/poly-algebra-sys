# Makefile for minimal polynomial project

CC = gcc
CFLAGS = -std=gnu99
LOADLIBES = -lm
SRC = minpoly.c subset_sum.c polynomial.c matrices.c integers.c roots.c interpolate.c resultant.c factoring.c algebraics.c calc_interface.c calculator.c
OBJ = minpoly.o subset_sum.o polynomial.o matrices.o integers.o roots.o interpolate.o resultant.o factoring.o algebraics.o calc_interface.o calculator.o
EXEC = minpoly subset_sum polynomial matrices integers roots interpolate resultant factoring algebraics calc_interface calculator

calculator: ${OBJ}

calc_interface: minpoly.o subset_sum.o polynomial.o matrices.o integers.o roots.o interpolate.o resultant.o factoring.o algebraics.o calc_interface.o

minpoly: CFLAGS += -Wall -DTEST_MINPOLY
minpoly: minpoly.o subset_sum.o polynomial.o roots.o integers.o

subset_sum: CFLAGS += -Wall -DTEST_SUBSET_SUM
subset_sum: subset_sum.o

polynomial: CFLAGS += -Wall -DTEST_POLYNOMIAL
polynomial: polynomial.o integers.o

matrices: CFLAGS += -Wall -DTEST_MATRICES
matrices: matrices.o

integers: CFLAGS += -Wall -DTEST_INTEGERS
integers: integers.o

roots: CFLAGS += -Wall -DTEST_ROOTS
roots: roots.o polynomial.o integers.o

interpolate: CFLAGS += -Wall -DTEST_INTERPOLATE
interpolate: interpolate.o polynomial.o matrices.o integers.o

resultant: CFLAGS += -Wall -DTEST_RESULTANT
resultant: resultant.o polynomial.o matrices.o integers.o interpolate.o

factoring: factoring.o

algebraics: CFLAGS += -Wall -DTEST_ALGEBRAICS
algebraics: algebraics.o minpoly.o subset_sum.o polynomial.o matrices.o integers.o roots.o interpolate.o resultant.o factoring.o

# remove object files prior to compiling test versions
test:
	rm -f ${OBJ}

depend:
	./remove_depends.sh
	gcc -MM ${SRC} >> Makefile

clean:
	rm -f ${OBJ} ${EXEC}

# dependencies listed by gcc -MM
minpoly.o: minpoly.c minpoly.h polynomial.h precision.h roots.h \
 subset_sum.h
subset_sum.o: subset_sum.c subset_sum.h precision.h
polynomial.o: polynomial.c polynomial.h precision.h integers.h
matrices.o: matrices.c matrices.h precision.h
integers.o: integers.c integers.h precision.h
roots.o: roots.c roots.h polynomial.h precision.h
interpolate.o: interpolate.c interpolate.h polynomial.h precision.h \
 matrices.h integers.h
resultant.o: resultant.c resultant.h polynomial.h precision.h matrices.h \
 interpolate.h
factoring.o: factoring.c factoring.h polynomial.h precision.h roots.h
algebraics.o: algebraics.c algebraics.h roots.h polynomial.h precision.h \
 minpoly.h resultant.h factoring.h
calc_interface.o: calc_interface.c calc_interface.h algebraics.h roots.h \
 polynomial.h precision.h
calculator.o: calculator.c calc_interface.h algebraics.h roots.h \
 polynomial.h precision.h
