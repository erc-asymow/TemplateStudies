GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O3 -I/usr/include/boost/ -L/usr/lib64/ -lboost_program_options

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
