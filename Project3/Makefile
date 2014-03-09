################################################################
#  Makefile for project 3
#
#  Will Spurgin
#  High Performance Scientific Computing
#  Math 3316
#  8 October 2013
################################################################

# compiler & flags
CXX = g++
CXXFLAGS = -O2 -std=c++0x

# makefile targets
all : test_lagrange.out test_newton.out compare.out test_lagrange2D.out runge_regular.out runge_chebyshev.out

test_lagrange.out : lagrange.cpp test_lagrange.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_newton.out : test_newton.cpp newton_interp.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

compare.out : compare.cpp newton_interp.cpp lagrange.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_lagrange2D.out : test_lagrange2D.cpp lagrange2D.cpp lagrange.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

runge_regular.out : runge_regular.cpp lagrange2D.cpp lagrange.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

runge_chebyshev.out : runge_chebyshev.cpp lagrange2D.cpp lagrange.cpp mat.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.out *~


####### End of Makefile #######
