###################################################################
#  Makefile for project 4
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 3316
#  6 November 2013
###################################################################

# compiler & flags
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	CXX = clang++
else
	CXX = g++
endif
CXXFLAGS = -O2


# makefile targets
all : test_Gauss2.out test_int.out test_adapt.out test_carbon.out application.out

test_Gauss2.out : test_Gauss2.cpp composite_Gauss2.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_int.out : test_int.cpp composite_int.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_adapt.out : test_adapt.cpp adaptive_int.cpp composite_int.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@	

test_carbon.out : test_carbon.cpp mat.cpp adaptive_int.cpp composite_int.cpp carbon.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

application.out : application.cpp mat.cpp adaptive_int.cpp composite_int.cpp carbon.cpp fd_newton.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.out *~


####### End of Makefile #######
