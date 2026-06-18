CXX ?= g++
CXXFLAGS ?= -std=c++17 -O1 -g

galois.o:
	$(CXX) $(CXXFLAGS) -c -o $@ galois.cpp
