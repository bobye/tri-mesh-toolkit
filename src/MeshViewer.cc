#include "MeshViewer.h"

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

MeshPainter::~MeshPainter() {
  delete [] vertex_array;
  delete [] normal_array;
  delete [] index_array;
}

MeshPainter::MeshPainter(TriMesh *pmesh) : obj(pmesh){

  coord_min_x = pmesh->coord_min_x;
  coord_min_y = pmesh->coord_min_y;
  coord_min_z = pmesh->coord_min_z;
  coord_max_x = pmesh->coord_max_x;
  coord_max_y = pmesh->coord_max_y;
  coord_max_z = pmesh->coord_max_z;

  vn = pmesh->vertex_num;
  fn = pmesh->facet_num;

  vertex_array = new float[3*vn];
  normal_array = new float[3*vn];
  index_array = new int[3*fn];
  
  for (int i=0; i<vn; i++) {
    Point p=pmesh->IV[i]->point();
    Vector n=pmesh->vertex_norm[i];
    vertex_array[3*i] = p.x();
    vertex_array[3*i+1] = p.y();
    vertex_array[3*i+2] = p.z();
    normal_array[3*i] = n.x();
    normal_array[3*i+1] = n.y();
    normal_array[3*i+2] = n.z();

  }

  for (int i=0; i<fn; i++) {
    HF_circulator hc = pmesh->IF[i]->facet_begin();
    int j=0;
    do{
      index_array[3*i+j] = hc->vertex()->index;
      ++j;
    }while(++hc != pmesh->IF[i]->facet_begin());        
  }

  LIST_NAME =glGenLists(1);


}

MeshPainter::MeshPainter(TriMesh *pmesh, Scalar_Fun *psfun) : obj(pmesh) {
}

void MeshPainter::draw(){  

  glColor3f(0.9, 0.9, 0.9);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glCullFace(GL_BACK);
  

  glDrawElements(GL_TRIANGLES, 3*fn, GL_UNSIGNED_INT, index_array);

  
}


void MeshPainter::prepare(){
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  //load geometric information from painte

  glNormalPointer(GL_FLOAT, 0, normal_array);
  glVertexPointer(3, GL_FLOAT, 0, vertex_array); 


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

MeshViewer::MeshViewer(int argc, char** argv) 
  :width(800), height(800) {
  currentMeshViewer = this; //static member need definition

  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowPosition(200.0, 0.0);
  glutInitWindowSize(width, height);
  glutCreateWindow("MeshTK - Viewer");

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
    

  glEnable(GL_LIGHTING);

  //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glEnable(GL_DEPTH_TEST);



}

void MeshViewer::display(){
  int n = currentMeshViewer->Painters.size();
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers

  for (int i=0; i<n; i++) 
    { 

      glPushMatrix();//push i-th matrix
      glCallList(currentMeshViewer->Painters[i]->LIST_NAME);

      
      //currentMeshViewer->Painters[i]->prepare();    
      //currentMeshViewer->Painters[i]->draw();
      glPushMatrix();//pop i-th matrix

    }
  glutSwapBuffers();
      
}






void MeshViewer::init(){

  
  /* 3D configuration */
  // clear the window color and depth buffer
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClearDepth(1.0);

  double center_x = (coord_min_x + coord_max_x) /2.;
  double center_y = (coord_min_y + coord_max_y) /2.;
  double center_z = (coord_min_z + coord_max_z) /2.;
  double length_z = coord_max_z - coord_min_z;
  
  double radio_x = (coord_max_x - coord_min_x) / (double) width;
  double radio_y = (coord_max_y - coord_min_y) / (double) height;
  double radio = (radio_x > radio_y)? radio_x: radio_y;

  //std::cout<< length_z << std::endl;
  glMatrixMode(GL_PROJECTION);

  glOrtho( center_x - radio * width/perfect_factor, center_x + radio * width/perfect_factor,
	   center_y - radio * height/perfect_factor, center_y + radio * height/perfect_factor,
	   //-10, 10);
	   center_z -  2* length_z , center_z +  2* length_z);//(NEW) set up our viewing area


  //glOrtho(-1,1,-1,1,-10,10);
  glMatrixMode(GL_MODELVIEW);


  /////////////////////////
  add_lights();
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


  painter->prepare();    

  glNewList(painter->LIST_NAME, GL_COMPILE_AND_EXECUTE);
  painter->draw();
  glEndList();


}

void MeshViewer::view(){

  glutMainLoop();
}













