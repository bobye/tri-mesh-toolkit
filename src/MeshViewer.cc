#include "MeshViewer.h"
#include <GL/glut.h>

MeshPainter::MeshPainter(TriMesh *pmesh){

}

MeshPainter::MeshPainter(TriMesh *pmesh, Scalar_Fun *psfun){
}

void MeshPainter::draw(){
}

MeshMarker::MeshMarker(TriMesh *pmesh, Bool_Fun* pbfun)
  :MeshPainter(pmesh) {
  
}

MeshMarker::MeshMarker(TriMesh *pmesh, Scalar_Fun *psfun, Bool_Fun *pbfun)
  :MeshPainter(pmesh, psfun) {
}


void MeshMarker::draw(){
}





MeshViewer::MeshViewer(int w, int h) 
  :width(w), height(h) {
  MeshViewer* currentMeshViewer = this; //static member need definition
}

void MeshViewer::display(){
  int n = currentMeshViewer->Painters.size();
  for (int i=0; i<n; i++) currentMeshViewer->Painters[i]->draw();
}




void MeshViewer::init(int argc, char** argv){

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition(200.0, 0.0);
  glutInitWindowSize(width, height);
  glutCreateWindow("MeshTK - Viewer");
  
  /* 3D configuration */
  /*
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearDepth(1.0);
    
  //glOrtho(-2,2,-2,2,-50,50);//(NEW) set up our viewing area
  
  glShadeModel(GL_SMOOTH);// Enable Smooth Shading
    
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
  //glLineWidth(1.0);
  //
  //glEnable(GL_CULL_FACE);
  //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
    
  //glutIdleFunc(idle);
    
  //lightsource(); light source configuration
  glEnable(GL_LIGHTING);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  //glEnable(GL_LIGHT0);//lighting
  //glEnable(GL_LIGHT1);
  //glEnable(GL_LIGHT2);
  glEnable(GL_DEPTH_TEST);
  */
  /////////////////////////

  glutDisplayFunc(display);
  


}

void MeshViewer::add_painter(MeshPainter *painter){
  Painters.push_back(painter);
}

void MeshViewer::view(){
  
}


