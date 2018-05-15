CXXFLAGS = -std=c++11
MPICXX=mpicxx

example: main_sequential.cpp
	${MPICXX} $(CXXFLAGS) -o sequential main_sequential.cpp physics.cpp

allclean:
	-rm *.o sequential
