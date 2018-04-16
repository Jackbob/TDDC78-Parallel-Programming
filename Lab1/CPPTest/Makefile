
EXECS=blurc thresc
MPICXX?=mpicxx
CXXFLAGS= -std=c++17

all: ${EXECS}

blurc: blurmain.cpp
	${MPICXX} ${CXXFLAGS} -o blurc blurmain.cpp ppmio.cpp gaussw.cpp blurfilter.cpp

thresc: thresmain.cpp
	${MPICXX} ${CXXFLAGS} -o thresc thresmain.cpp thresfilter.cpp ppmio.cpp

clean:
	rm -f ${EXECS}