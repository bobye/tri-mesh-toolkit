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

LIBDIR=lib
LIBPATH = -L$(LIBDIR)

OBJDIR =obj

OBJECTS = TriMesh_Base.o TriMesh_Curv.o TriMesh_SIFT.o TriMesh_UI.o DynamicTriMesh.o MeshViewer.o mesh_assist.o ##temp objects

VPATH = src src/meshtk



default: all

$(OBJDIR)/%.o: %.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $<

lib: $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(CPP) -shared -Wl,-soname,$(LIBDIR)/libmeshtk.so -o $(LIBDIR)/libmeshtk.so $^

MeshTK: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) $(LIBTCLAP) -lmeshtk -o MeshTK $< 

test: test.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -lmeshtk -o test $<

all: lib MeshTK test

clean: 
	$(RM) $(OBJDIR)/*.o $(LIBDIR)/*.so MeshTK




