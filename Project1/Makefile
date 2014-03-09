#####################################
#  Will Spurgin 
#  Math 3316
#  14 September 2013
#####################################

CXX = g++
CXXFLAGS = -O0 -g

all : taylor_series.out floating_error.out

floating_error.out : floating_error.cpp mat.o
	$(CXX) $(CXXFLAGS) $^  -o $@

taylor_series.out : taylor_series.cpp mat.o
	$(CXX) $(CXXFLAGS) $^  -o $@

mat.o : mat.cpp mat.h
	$(CXX) $(CXXFLAGS) -c mat.cpp -o $@

clean :
	\rm -f *.o *.out

veryclean : clean
	\rm -f *.txt *~


##### End of Makefile #####
