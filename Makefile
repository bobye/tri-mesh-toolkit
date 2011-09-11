CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include -lCGAL
VPATH = src include
BIN = bin


default: MeshTK

%.o: %.cc
	$(CPP) -c $(CPPFLAGS) -o $@ $<

MeshTK: main.o
	$(CPP) $(CPPFLAGS) -o $(BIN)/MeshTK $^

clean: main.o
	$(RM) $^
