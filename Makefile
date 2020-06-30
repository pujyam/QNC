#CXX = c++
CXX = clang++ 
CXXFLAGS+=-stdlib=libc++ -nostdinc++
#CXX = g++
#CXX = /opt/local/bin/g++-mp-6

INC_DIR:=.

ifeq ($(OS),Windows_NT)
    #CCFLAGS +=-D WIN32
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        #LIB_FLAGS+=-lopenblas
        #EXTRA_OPT=-fwhole-program
    endif
    ifeq ($(UNAME_S),Darwin)
        #LIB_FLAGS+=-framework Accelerate
    endif
endif

FLAGS=-std=c++17 -I$(INC_DIR)/include -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK

## As the Armadillo library uses recursive templates, compilation times depend on level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups, but compilation takes longer
## -O3: enables auto vectorisation when using gcc
#OPT=-O0 -g -Wall -pedantic
OPT=-O3 -Wall -pedantic

## Uncomment line if you're compiling all source files into one program in a single hit

#DEBUG=-DVERY_VERBOSE
#DEBUG+=-DARMA_EXTRA_DEBUG
## Uncomment the above line(s) to enable low-level debugging.
## Lots of debugging information will be printed when a compiled program is run.
## Please enable this option when reporting bugs.

FINAL=-DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## Not recommended unless your code has been first thoroughly tested!

CXXFLAGS+=$(OPT) $(EXTRA_OPT) $(FLAGS) $(DEBUG) $(FINAL)

all: qc

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

qc: qgate.o qunitary.o qperm.o qasm.o qc.o
	$(CXX) -o $@ $^ $(LIB_FLAGS)

.PHONY: clean
clean:
	rm -rf *.o qc qc.dSYM
