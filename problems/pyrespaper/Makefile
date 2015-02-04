# This Makefile compiles the shared dynamic library to access the IAS15 integrator.
CC=cc 
OPT=-std=c99 -O3 -DCOLLISIONS_NONE

all: 
	$(CC) $(OPT) -c -fPIC -Wall ../../src/particle.c
	$(CC) $(OPT) -c -fPIC -I../../src/ -Wall problem_tools.c
	$(CC) $(OPT) -shared -o problem_tools.so particle.o problem_tools.o

	@echo ""	
	@echo "The shared library problem_tools.so has been created. Cleaning up"
	@echo ""	

	rm -f *.o

