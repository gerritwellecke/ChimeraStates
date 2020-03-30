 ###############################################################################
 # Author:   Gerrit Wellecke
 # Project:  Chimera states in populations of Kuramoto oscillators
 # Date:     March 2020
 # 
 # Generic makefile for this project. This will work if executed on a Darwin
 # system. The flags for compiling openMP need to be changed for Linux systems.
 ###############################################################################

BINARY := kuramoto_simulation.out
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
HEADER := $(wildcard *.hpp)

CXX := g++
CXXFLAGS := -Wall -std=c++2a -O2 -Xpreprocessor -fopenmp
LFLAGS := -lm -lomp


# link binary from object files
$(BINARY): $(OBJECTS)
	$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -o $@

# compile all .cpp files in directory
%.o: %.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) $< -c

# remove all GCC compiled files, the main binary and all written data
.PHONY: clean
clean: 
	@echo "removing all object files"
	@rm -f *.o
	@echo "removing all data files"
	@rm -f *.dat
