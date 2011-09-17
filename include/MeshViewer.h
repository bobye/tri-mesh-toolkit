#ifndef _MESHVIEWER_H_
#define _MESHVIEWER_H_

#include "TriMesh.h"

class MeshPainter {
  TriMesh *obj;
  Scalar_Fun *color;

  

public:
  MeshPainter(TriMesh *, Scalar_Fun *); // with color ramping 
  MeshPainter(TriMesh *); 
  void draw();
}

class MeshMarker : public MeshPainter {
  
  Bool_Fun *mark;

public:
  MeshMarker(TriMesh *, Scalar_Fun *, Bool_Fun *);
  MeshMarker(TriMesh *, Bool_Fun *);
  void draw();
}

class MeshViewer {
  vector<MeshPainter*> Painters;

public:
  MeshViewer(int ,int );//init windows and gl setting
  void add_painter(MeshPainter *);
  
  void view();//main loop
}

#endif /* _MESHVIEWER_H_ */
