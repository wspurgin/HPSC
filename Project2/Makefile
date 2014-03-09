#########################################################
#  Makefile
#
#  Will Spurgin 
#  High Performance Scientific Computing 
#  Math 3316
#  23 September 2013
#########################################################

all : application.out test_bisect.out test_newton.out test_fd_newton.out

application.out : application.cpp fd_newton.cpp newton.cpp bisection.cpp
	g++ $^ -o $@

test_fd_newton.out : test_fd_newton.cpp fd_newton.cpp
	g++ $^ -o $@

test_newton.out : test_newton.cpp newton.cpp
	g++ $^ -o $@

test_bisect.out : test_bisect.cpp bisection.cpp
	g++ $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.out *~


##### End of Makefile #####
