## Please modify the following options 

CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include 

LIBOPT = -lCGAL -lglut -lGL -lGLU -I/usr/local/include/tclap/


## DO NOT modify the following lines.

VPATH = src include



default: MeshTK

%.o: %.cc
	$(CPP) -c $(CPPFLAGS) $(LIBOPT) -o $@ $<

MeshTK: main.o  TriMesh.o MeshViewer.o mesh_assist.o
	$(CPP) $(CPPFLAGS) $(LIBOPT) -o MeshTK $^

clean: 
	$(RM) *.o MeshTK
