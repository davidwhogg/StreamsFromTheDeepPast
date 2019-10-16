RM= /bin/rm -vf

all: EM

EM: EM.o
	gcc -o EM EM.o -lm -lgsl -lgslcblas

EM.o: EM.c
	gcc -Wall -c EM.c

EM.so: EM.c
	gcc -fpic -Wall -c EM.c
	gcc -shared -o EM.so EM.o -lm -lgsl -lgslcblas

.PHONY: clean

clean:
	$(RM) EM.o EM.so
