GCC=g++
#GCC=/cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0/x86_64-centos7/bin/g++
CXXFLAGS=`root-config --libs --cflags` -O3 -I/usr/include/boost/ -L/usr/lib64/ -lboost_program_options -I/home/users/bianchini/Eigen/eigen-3.4.0

SRCDIR=.
BINDIR=.

OBJ=main.cpp

.PHONY: all
all:
	$(info, "--- Full compilation --- ")	
	$(info, "-> if you want just to recompile something use 'make fast' ")	
	$(info, "------------------------ ")	
	$(MAKE) main.o

.PHONY: fast
fast:
	$(MAKE) main.o

main.o: main.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/main $(OBJ)

jac.o: jac.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/jac jac.cpp  

fit.o: fit.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/fit fit.cpp  
