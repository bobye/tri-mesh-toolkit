/*
  FILE: meshtk/MeshViewer.hh This file is part of MeshTK.
  It is a C++ header file which defines the interface of class 
  MeshPainter, MeshRamper, MeshMarker, and MeshViewer. 
  
  Copyright (C) 2011 Jianbo YE

  MeshTK is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  MeshTK is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA  
*/

#ifndef _MESHVIEWER_HH_
#define _MESHVIEWER_HH_

#include <vector>
#include "TriMesh.hh"
#include <GL/glut.h>

namespace meshtk {

  class MeshPainter {

  protected:
    TriMesh *obj;
    GLuint vn; //vertex number
    GLuint fn; //facet number
  


    GLfloat *vertex_array;
    GLfloat *normal_array;
    GLuint *index_array;

  public:
    GLuint LIST_NAME; //name of display list

    GLfloat coordinate_min_x, coordinate_min_y, coordinate_min_z;
    GLfloat coordinate_max_x, coordinate_max_y, coordinate_max_z;


    MeshPainter(TriMesh *); 
    virtual void prepare();
    virtual void draw();
  
    ~MeshPainter();
    
  };

  class MeshRamper : public MeshPainter {

    GLfloat *color_array;
  public:
    MeshRamper(TriMesh *, ScalarFunction *);
    void prepare();

    ~MeshRamper();
  };

  class MeshMarker : public MeshPainter {
    GLuint mark_count;
    GLuint *mark_array;

  public:
    MeshMarker(TriMesh *, BooleanFunction *);
    void prepare();
    void draw();
    ~MeshMarker();
  };

  class MeshViewer {
    int width, height;
    GLfloat coordinate_min_x, coordinate_min_y, coordinate_min_z;
    GLfloat coordinate_max_x, coordinate_max_y, coordinate_max_z;

    GLfloat center_x, center_y, center_z, length_z;
    std::vector<MeshPainter*> Painters;
  
    static MeshViewer* currentMeshViewer;
    static void display();
    static void reshape(int, int);
    static void keyboard(unsigned char, int, int);
    static void motion(int, int);
    static void mouse(int, int, int, int);




    void add_lights();  

  public:
    MeshViewer(int ,char** );//init windows and gl setting
    void init(); //call in main()

    void add_painter(MeshPainter *);
    void view();//main loop

  };
}
#endif /* _MESHVIEWER_HH_ */
