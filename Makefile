CXX=g++
CC=$(CXX)

#
CFLAGS= -O2 -march=native -pipe -funroll-loops
CXXFLAGS=$(CFLAGS)

OBJS= mainRDME.o
EXE= rdme

all:mainRDME
	mv mainRDME $(EXE)
main:$(OBJS)

clean:
	rm -f $(OBJS) $(EXE)

redo: clean all
