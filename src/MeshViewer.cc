#include "MeshViewer.h"
#include <GL/glut.h>
#include "TriMesh.h"
#include "mesh_topo.h"

MeshPainter::MeshPainter(TriMesh *pmesh) : obj(pmesh){

}

MeshPainter::MeshPainter(TriMesh *pmesh, Scalar_Fun *psfun) : obj(pmesh) {
}

void MeshPainter::draw(){  
  
  glColor3f(0.9, 0.9, 0.9);
  
  for(Facet_iterator fi = obj->P.facets_begin(); fi != obj->P.facets_end(); ++fi){
    HF_circulator hc = fi->facet_begin();
        
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_LINE_LOOP);    
    
    do{
      glVertex3f(hc->vertex()->point().x(), hc->vertex()->point().y(), hc->vertex()->point().z());
    }while(++hc != fi->facet_begin());
    glEnd();
        
  }

}

MeshMarker::MeshMarker(TriMesh *pmesh, Bool_Fun* pbfun)
  :MeshPainter(pmesh) {
  
}

MeshMarker::MeshMarker(TriMesh *pmesh, Scalar_Fun *psfun, Bool_Fun *pbfun)
  :MeshPainter(pmesh, psfun) {
}


void MeshMarker::draw() {
  MeshPainter::draw();
  
}



MeshViewer* MeshViewer::currentMeshViewer;

MeshViewer::MeshViewer(int w, int h) 
  :width(w), height(h) {
  currentMeshViewer = this; //static member need definition
}

void MeshViewer::display(){
  int n = currentMeshViewer->Painters.size();
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers

  for (int i=0; i<n; i++) 
    { 
      glPushMatrix();//push i-th matrix
      currentMeshViewer->Painters[i]->draw();
      glPushMatrix();//pop i-th matrix
    }
  glutSwapBuffers();
      
}






void MeshViewer::init(int argc, char** argv){

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowPosition(200.0, 0.0);
  glutInitWindowSize(width, height);
  glutCreateWindow("MeshTK - Viewer");
  
  /* 3D configuration */

  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClearDepth(1.0);

  glOrtho(-1,1,-1,1,-10,10);//(NEW) set up our viewing area
  
  glShadeModel(GL_SMOOTH);// Enable Smooth Shading
    
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
  //glLineWidth(1.0);
  //
  glEnable(GL_CULL_FACE);
  //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
    
  //glutIdleFunc(idle);
    
  //lightsource(); light source configuration
  //glEnable(GL_LIGHTING);
  //glCullFace(GL_BACK);
  //glEnable(GL_CULL_FACE);

  //glEnable(GL_LIGHT0);//lighting
  //glEnable(GL_LIGHT1);
  //glEnable(GL_LIGHT2);
  glEnable(GL_DEPTH_TEST);

  /////////////////////////

  glutDisplayFunc(display);




}

void MeshViewer::add_painter(MeshPainter *painter){
  Painters.push_back(painter);
}

void MeshViewer::view(){
  glutMainLoop();
}


