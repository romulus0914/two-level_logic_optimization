GCC = gcc
GXX = g++

GCCFLAG = -std=c99 -O3
GXXFLAG = -std=c++11 -O3
LIB = -lm -fopenmp

CSOURCE = two-level_logic_opt.c
CPPSOURCE = two-level_logic_opt.cpp
CEXE = two-level_logic_opt
CPPEXE = two-level_logic_opt_cpp

all: c cpp

c:
	$(GCC) $(GCCFLAG) $(LIB) $(CSOURCE) -o $(CEXE)

cpp:
	$(GXX) $(GXXFLAG) $(LIB) $(CPPSOURCE) -o $(CPPEXE)

clean:
	rm  -rf $(CEXE) $(CPPEXE)
