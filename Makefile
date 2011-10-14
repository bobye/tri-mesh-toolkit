## Please modify the following options properly 

## make sure your environment variable PETSC_DIR are correctly given
PETSC_DIR := /home/bobye/pub/petsc/petsc-3.2-p3
include ${PETSC_DIR}/conf/variables

## options for CGAL, OpenGL and freeglut 
LIBOPT = -lCGAL -lglut -lGL -lGLU 
## option for TCLAP 
LIBTCLAP = -I/usr/local/include/tclap 

## add this directory to configuration file of ld
## it may be /etc/ld.so.conf, and run ldconfig -v to update cache as root.
MESHTKLIBPATH = /opt/lib

## executable toolkit directory, where you do the benchworks
MESHTKPATH = ~/data/meshtk_workshop

CC = gcc
CPP = g++
AR = ar


#####################################################################################
## DO NOT change anything below here
CPPFLAGS = -Wall -I include

TOPDIR = $(shell pwd)

## temporary directory to store dynamic libraries
LIBDIR= lib
LIBPATH = -L$(LIBDIR)

## temporary directory to store object files
OBJDIR = obj

## objects to build the library
OBJECTS = TriMesh_Base.o TriMesh_Curv.o TriMesh_SIFT.o TriMesh_Dist.o TriMesh_UI.o TriMesh_PETSc.o DynamicTriMesh.o MeshViewer.o mesh_assist.o ##temp objects

## directory of example executable
EXDIR = example

## examples
EXAMPLE = ex0

## c++ source files
VPATH = src src/meshtk src/example


## rules start 

default: lib

## library
$(OBJDIR)/%.o: %.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $< 

$(OBJDIR)/TriMesh_PETSc.o: src/meshtk/TriMesh_PETSc.cc
	$(CPP) -c -fPIC $(CPPFLAGS) $(LIBOPT) -o $@ $< ${PETSC_CCPPFLAGS} 


lib: $(addprefix $(OBJDIR)/, $(OBJECTS)) 
	$(CPP) -shared -Wl,-soname,libmeshtk.so.1 -o $(LIBDIR)/libmeshtk.so.1.0 $^ ${PETSC_MAT_LIB}
	cp $(LIBDIR)/*.so* $(MESHTKLIBPATH)
	ln -sf $(MESHTKLIBPATH)/libmeshtk.so.1.0 $(MESHTKLIBPATH)/libmeshtk.so.1 
	ln -sf $(MESHTKLIBPATH)/libmeshtk.so.1.0 $(MESHTKLIBPATH)/libmeshtk.so 

## toolkits
MeshTK: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) $(LIBTCLAP) -L$(MESHTKLIBPATH) -lmeshtk -o $(MESHTKPATH)/MeshTK $<
test: main.cc
	$(CPP) $(CPPFLAGS) $(LIBOPT) -L$(MESHTKLIBPATH) -lmeshtk -o $(MESHTKPATH)/test $<


toolkit: MeshTK test

## examples
$(EXDIR)/ex0: ex0.cc 
	$(CPP) $(CPPFLAGS) $(LIBOPT) -L$(MESHTKLIBPATH) -lmeshtk -o $@ $< 

example: $(addprefix $(EXDIR)/, $(EXAMPLE))

## 
## debug
all: lib example toolkit


clean:	
	$(RM) $(OBJDIR)/*.o $(LIBDIR)/*.so* $(LIBDIRS)/*.a $(EXDIR)/ex* $(MESHTKPATH)/MeshTK $(MESHTKPATH)/test














