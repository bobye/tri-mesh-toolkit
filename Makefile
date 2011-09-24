## Please modify the following options 

CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include 

LIBOPT = -lCGAL -lglut -lGL -lGLU 

## DO NOT change anything below here
LIBDIR=lib
LIBPATH = -L$(LIBDIR)

OBJDIR =obj

OBJECTS = TriMesh.o DynamicTriMesh.o MeshViewer.o mesh_assist.o ##temp objects

VPATH = src src/meshtk



default: all

$(OBJDIR)/%.o: %.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $<

lib: $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(CPP) -shared -Wl,-soname,$(LIBDIR)/libmeshtk.so -o $(LIBDIR)/libmeshtk.so $^

MeshTK: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -I/usr/local/include/tclap/ -lmeshtk -o MeshTK $< 

all: lib MeshTK

clean: 
	$(RM) $(OBJDIR)/*.o $(LIBDIR)/*.so MeshTK

