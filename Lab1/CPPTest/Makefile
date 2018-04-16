
CFLAGS = 

LFLAGS= -lpthread -lm

all: blurc thresc

clean:
	-$(RM) *.o blur thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	mpicxx -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	mpicxx -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cf - *.cpp *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
