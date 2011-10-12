## Please modify the following options properly 

## make sure your environment variable PETSC_DIR are correctly given
PETSC_DIR := /home/bobye/pub/petsc/petsc-3.2-p3
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

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

OBJECTS = TriMesh_Base.o TriMesh_Curv.o TriMesh_SIFT.o TriMesh_Dist.o TriMesh_UI.o TriMesh_PETSc.o DynamicTriMesh.o MeshViewer.o mesh_assist.o ##temp objects

EXDIR = example

EXAMPLE = ex0

VPATH = src src/meshtk src/example


## library
$(OBJDIR)/%.o: %.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $<

$(OBJDIR)/TriMesh_PETSc.o: TriMesh_PETSc.cc
	-${CLINKER} -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $(OBJDIR)/TriMesh_PETSc.o src/meshtk/TriMesh_PETSc.cc ${PETSC_CCPPFLAGS}

libmeshtk: $(addprefix $(OBJDIR)/, $(OBJECTS)) 
	-${CLINKER} -shared -Wl,-soname,$(LIBDIR)/libmeshtk.so -o $(LIBDIR)/libmeshtk.so $^ ${PETSC_MAT_LIB}

## toolkits
MeshTK: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) $(LIBTCLAP) -lmeshtk -o MeshTK $< 

test: test.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -lmeshtk -o test $< 

## examples
$(EXDIR)/ex0: ex0.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBPATH) -lmeshtk -o $@ $< 

example: $(addprefix $(EXDIR)/, $(EXAMPLE))

all: libmeshtk MeshTK test example

cleanall:	
	$(RM) $(OBJDIR)/*.o $(LIBDIR)/*.so $(EXDIR)/ex* MeshTK test














