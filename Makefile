## Please modify the following options properly 

CC = gcc
CPP = g++
CPPFLAGS = -Wall -I include

## options for CGAL, OpenGL and freeglut 
LIBOPT = -lCGAL -lglut -lGL -lGLU 
## option for TCLAP 
LIBTCLAP = -I/usr/local/include/tclap

#####################################################################################
## DO NOT change anything below here

LIBDIR= lib
LIBPATH = -L$(LIBDIR)

OBJDIR = obj

OBJECTS = TriMesh_Base.o TriMesh_Curv.o TriMesh_SIFT.o TriMesh_Dist.o TriMesh_UI.o DynamicTriMesh.o MeshViewer.o mesh_assist.o ##temp objects

EXDIR = example

EXAMPLE = ex0

VPATH = src src/meshtk src/example



default: all
## library
$(OBJDIR)/%.o: %.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $<

lib: $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(CPP) -shared -Wl,-soname,$(LIBDIR)/libmeshtk.so -o $(LIBDIR)/libmeshtk.so $^

## toolkits
MeshTK: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) $(LIBTCLAP) -lmeshtk -o MeshTK $< 

test: test.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -lmeshtk -o test $<

## examples
$(EXDIR)/ex0: ex0.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -lmeshtk -o $@ $<

example: $(addprefix $(EXDIR)/, $(EXAMPLE))

all: lib MeshTK test

clean: 
	$(RM) $(OBJDIR)/*.o $(LIBDIR)/*.so $(EXDIR)/* MeshTK test














