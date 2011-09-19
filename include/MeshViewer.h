#ifndef _MESHVIEWER_H_
#define _MESHVIEWER_H_

#include <vector>
#include "TriMesh.h"


class MeshPainter {
  TriMesh *obj;
  Scalar_Fun *color;
  //  double mat_ambient[4], mat_diffuse[4], mat_specular[4];

public:
  double coord_min_x, coord_min_y, coord_min_z;
  double coord_max_x, coord_max_y, coord_max_z;

  MeshPainter(TriMesh *, Scalar_Fun *); // with color ramping 
  MeshPainter(TriMesh *); 
  void draw();

};

class MeshMarker : public MeshPainter {
  
  Bool_Fun *mark;

public:
  MeshMarker(TriMesh *, Scalar_Fun *, Bool_Fun *);
  MeshMarker(TriMesh *, Bool_Fun *);
  void draw();
};

class MeshViewer {
  int width, height;
  double coord_min_x, coord_min_y, coord_min_z;
  double coord_max_x, coord_max_y, coord_max_z;

  std::vector<MeshPainter*> Painters;
  
  static MeshViewer* currentMeshViewer;
  static void display() ;
  
  void add_lights();  

public:
  MeshViewer(int ,int );//init windows and gl setting
  void init(int, char**); //call in main()

  void add_painter(MeshPainter *);
  void view();//main loop



};

#endif /* _MESHVIEWER_H_ */
