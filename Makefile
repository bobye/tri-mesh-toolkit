## Please modify the following options 

CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include 

LIBOPT = -lCGAL  -I/usr/local/include/tclap/


## DO NOT modify the following lines.

VPATH = src include



default: MeshTK

%.o: %.cc
	$(CPP) -c $(CPPFLAGS) $(LIBOPT) -o $@ $<

MeshTK: main.o mesh_update.o mesh_vis.o
	$(CPP) $(CPPFLAGS) $(LIBOPT) -o MeshTK $^

clean: main.o
	$(RM) $^ MeshTK
