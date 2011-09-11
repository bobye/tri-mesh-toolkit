CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include -lCGAL
VPATH = src include



default: MeshTK

%.o: %.cc
	$(CPP) -c $(CPPFLAGS) -o $@ $<

MeshTK: main.o
	$(CPP) $(CPPFLAGS) -o MeshTK $^

clean: main.o
	$(RM) $^
