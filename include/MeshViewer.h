#ifndef _MESHVIEWER_H_
#define _MESHVIEWER_H_

#include <vector>
#include "TriMesh.h"
#include <GL/glut.h>


class MeshPainter {
  TriMesh *obj;
  GLuint vn; //vertex number
  GLuint fn; //facet number
  


  GLfloat *vertex_array;
  GLfloat *normal_array;
  GLuint *index_array;

public:
  GLuint LIST_NAME; //name of display list

  GLfloat coord_min_x, coord_min_y, coord_min_z;
  GLfloat coord_max_x, coord_max_y, coord_max_z;


  MeshPainter(TriMesh *); 
  virtual void prepare();
  void draw();
  
  ~MeshPainter();
};

class MeshRamper : public MeshPainter {

  GLfloat *color;
public:
  MeshRamper(TriMesh *, Scalar_Fun *);
  void prepare();
};

class MeshMarker : public MeshPainter {
  
  GLuint *mark_array;

public:
  MeshMarker(TriMesh *, Bool_Fun *);
  void prepare();
};

class MeshViewer {
  int width, height;
  GLfloat coord_min_x, coord_min_y, coord_min_z;
  GLfloat coord_max_x, coord_max_y, coord_max_z;

  std::vector<MeshPainter*> Painters;
  
  static MeshViewer* currentMeshViewer;
  static void display() ;
  
  void add_lights();  

public:
  MeshViewer(int ,char** );//init windows and gl setting
  void init(); //call in main()

  void add_painter(MeshPainter *);
  void view();//main loop



};

#endif /* _MESHVIEWER_H_ */
