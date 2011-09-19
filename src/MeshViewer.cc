#include "MeshViewer.h"
#include <GL/glut.h>
#include "TriMesh.h"
#include "mesh_topo.h"


GLfloat light_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };



GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };


GLfloat mat_ambient[] = { .3, .5, .6, 1.0 };
GLfloat mat_diffuse[] = { .3, .5, .6, 1.0 };
GLfloat mat_specular[] = { .3, .5, .6, 1.0 };
GLfloat mat_shininess[] = {30};

const double perfect_factor = 1.414;


MeshPainter::MeshPainter(TriMesh *pmesh) : obj(pmesh){

  coord_min_x = pmesh->coord_min_x;
  coord_min_y = pmesh->coord_min_y;
  coord_min_z = pmesh->coord_min_z;
  coord_max_x = pmesh->coord_max_x;
  coord_max_y = pmesh->coord_max_y;
  coord_max_z = pmesh->coord_max_z;
  
}

MeshPainter::MeshPainter(TriMesh *pmesh, Scalar_Fun *psfun) : obj(pmesh) {
}

void MeshPainter::draw(){  
  
  glColor3f(0.9, 0.9, 0.9);
  glPolygonMode(GL_FRONT, GL_FILL);
  glCullFace(GL_BACK);

  
  glBegin(GL_TRIANGLES);    
  
  for(Facet_iterator fi = obj->P.facets_begin(); fi != obj->P.facets_end(); ++fi){
    HF_circulator hc = fi->facet_begin();
    
    do{
      Vertex_handle v = hc->vertex();
      Point p=v->point();
      Vector n=obj->vertex_norm[v->index];
      glNormal3f(n.x(), n.y(), n.z());
      glVertex3f(p.x(), p.y(), p.z());
    }while(++hc != fi->facet_begin());        
  }
  glEnd();


  //    glutSolidSphere(0.5, 40, 16);


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
  // clear the window color and depth buffer
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClearDepth(1.0);

  double center_x = (coord_min_x + coord_max_x) /2.;
  double center_y = (coord_min_y + coord_max_y) /2.;
  double center_z = (coord_min_z + coord_max_z) /2.;
  double length_z = center_z - coord_min_z;
  
  double radio_x = (coord_max_x - coord_min_x) / width;
  double radio_y = (coord_max_y - coord_min_y) / height;
  double radio = (radio_x > radio_y)? radio_x: radio_y;

  glOrtho( center_x - radio * width/perfect_factor, center_x + radio * width/perfect_factor,
	   center_y - radio * height/perfect_factor, center_y + radio * height/perfect_factor,
	   center_z - 2* length_z , center_z + 2* length_z);//(NEW) set up our viewing area


  //glOrtho(-1,1,-1,1,-10,10);

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
  glEnable(GL_LIGHTING);
  add_lights();

  //glEnable(GL_LIGHT0);//lighting
  //glEnable(GL_LIGHT1);
  //glEnable(GL_LIGHT2);
  glEnable(GL_DEPTH_TEST);

  /////////////////////////

  glutDisplayFunc(display);






}


void MeshViewer::add_lights(){

  light_position[0] = coord_max_x * perfect_factor;
  light_position[1] = coord_max_y * perfect_factor;
  light_position[2] = coord_max_z * perfect_factor;

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glEnable(GL_LIGHT0);

  glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);


}
void MeshViewer::add_painter(MeshPainter *painter){
  Painters.push_back(painter);

  if (painter->coord_min_x < coord_min_x) coord_min_x = painter->coord_min_x;
  if (painter->coord_min_y < coord_min_y) coord_min_y = painter->coord_min_y;
  if (painter->coord_min_z < coord_min_z) coord_min_z = painter->coord_min_z;
  if (painter->coord_max_x > coord_max_x) coord_max_x = painter->coord_max_x;
  if (painter->coord_max_y > coord_max_y) coord_max_y = painter->coord_max_y;
  if (painter->coord_max_z > coord_max_z) coord_max_z = painter->coord_max_z;

}

void MeshViewer::view(){

  glutMainLoop();
}




