
EXECS=blurc thresc
MPICXX?=mpicxx

all: ${EXECS}

blurc: blurmain.cpp
	${MPICXX} -o blurc blurmain.cpp ppmio.cpp gaussw.cpp blurfilter.cpp

thresc: thresmain.cpp
	${MPICXX} -o thresc thresmain.cpp thresfilter.cpp ppmio.cpp

clean:
	rm -f ${EXECS}