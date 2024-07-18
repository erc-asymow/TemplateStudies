GCC=g++
#GCC=/cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0/x86_64-centos7/bin/g++
#CXXFLAGS=`root-config --libs --cflags` -O3 -I/usr/include/boost/ -L/usr/lib64/ -lboost_program_options -I/home/users/bianchini/Eigen/eigen-3.4.0
CXXFLAGS=`root-config --libs --cflags` -lRooFit -lRooFitCore -lMinuit2 -O3 -I/usr/include/boost/ -L/usr/lib64/ -lboost_program_options -I/usr/include/eigen3

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

jac2: jac2.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/jac2 jac2.cpp  

jac3: jac3.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/jac3 jac3.cpp  

jac4: jac4.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/jac4 jac4.cpp  

fit: fit.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/fit fit.cpp  

fit_grid: fit_grid.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/fit_grid fit_grid.cpp  

debug.o: debug.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/debug debug.cpp  

estim: estim.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/estim estim.cpp

jac2tmp: jac2.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/jac2tmp jac2.cpp  

getparam: getparam.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/getparam getparam.cpp  

toy: toy.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/toy toy.cpp  

massfit: massfit.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/massfit massfit.cpp  

resolfit: resolfit.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/resolfit resolfit.cpp  

massscales: massscales.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/massscales massscales.cpp  

massscales_data: massscales_data.cpp
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/massscales_data massscales_data.cpp  
