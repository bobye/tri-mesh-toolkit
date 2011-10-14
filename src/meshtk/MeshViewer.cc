/*
  FILE: MeshViewer.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implements the classes related to visualization.
  
  Copyright (C) 2011 Jianbo YE

  tri-mesh-toolkit is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  tri-mesh-toolkit is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA  
*/


#include "meshtk/MeshViewer.hh"
#include "meshtk/mesh_assist.hh"

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_topo.hh"

namespace meshtk {


  MeshViewer* MeshViewer::currentMeshViewer;
  GLfloat scale = 1.;
  int origin_x, origin_y;
  int transform_x=0., transform_y=0.;
  //int displace_x=0, displace_y=0;
  bool left_button = false , right_button= false;
  GLfloat CTM[16];

  GLfloat light_ambient[] = { .4, .4, .4, 1.0 };
  GLfloat light_diffuse[] = { .8, .8, .8, 1.0 };
  GLfloat light_specular[] = { .5, .5, .5, 1.0 };



  GLfloat light_position0[] = { 1.0, 1.0, 1.0, 0.0 };
  GLfloat light_position1[] = { -1.0, -1.0, -1.0, .0};

  GLfloat mat_ambient[] = { .5, .5, .5, 1.0 };
  GLfloat mat_diffuse[] = { .5, .5, .5, 1.0 };
  GLfloat mat_specular[] = { .1, .1, .1, 1.0 };
  GLfloat mat_shininess[] = {50};

  const GLfloat perfect_factor = 1.414;

  MeshPainter::~MeshPainter() {
  }

  MeshPainter::MeshPainter(TriMesh *pmesh) : obj(pmesh){

    coordinate_min_x = pmesh->coordinate_min_x;
    coordinate_min_y = pmesh->coordinate_min_y;
    coordinate_min_z = pmesh->coordinate_min_z;
    coordinate_max_x = pmesh->coordinate_max_x;
    coordinate_max_y = pmesh->coordinate_max_y;
    coordinate_max_z = pmesh->coordinate_max_z;

    vn = pmesh->vertex_num;
    fn = pmesh->facet_num;
    vertex_array = & pmesh->vertex_array[0];
    normal_array = & pmesh->normal_array[0];
    index_array = & pmesh->tri_index_array[0];

    /*
    vertex_array = new GLfloat[3*vn];
    normal_array = new GLfloat[3*vn];
    index_array = new GLuint[3*fn];
  
    for (GLuint i=0; i<vn; i++) {
      Point p=pmesh->IV[i]->point();
      Vector n=pmesh->vertex_norm[i];
      vertex_array[3*i] = p.x();
      vertex_array[3*i+1] = p.y();
      vertex_array[3*i+2] = p.z();
      normal_array[3*i] = n.x();
      normal_array[3*i+1] = n.y();
      normal_array[3*i+2] = n.z();

    }

    for (GLuint i=0; i<fn; i++) {
      HF_circulator hc = pmesh->IF[i]->facet_begin();
      GLuint j=0;
      do{
	index_array[3*i+j] = hc->vertex()->index;
	++j;
      }while(++hc != pmesh->IF[i]->facet_begin());        
    }
    */
    LIST_NAME =glGenLists(1);


  }


  void MeshPainter::draw(){  
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.3, 0.5, 0.6);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glCullFace(GL_BACK);
  

    glDrawElements(GL_TRIANGLES, 3*fn, GL_UNSIGNED_INT, index_array);
    glDisable(GL_COLOR_MATERIAL);
  
  
  }


  void MeshPainter::prepare(){
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    //load geometric information from painte

    glNormalPointer(GL_DOUBLE, 0, normal_array);
    glVertexPointer(3, GL_DOUBLE, 0, vertex_array); 


  }





  void color_ramping(GLfloat *color, ScalarFunction* psfun, bool equalhist = false){
    ScalarFunction s = *psfun;
    GLuint size = s.size();

    if (equalhist) {
      std::vector<size_t> index(size);
      for (unsigned i=0; i<size; ++i) index[i]=i;

      std::sort(index.begin(), index.end(), index_cmp<std::vector<double>&>(s));
      double sum=0;
      for (unsigned i=0; i<size; ++i) 
	{
	  s[index[i]] = sum++;
	}
    }

    double vmax= *std::max_element(s.begin(), s.end());
    double vmin= *std::min_element(s.begin(), s.end());
    //double vmax= sum-1; double vmin =0;
    double dv = vmax - vmin;

    //  std::cout<<vmax <<"\t"<< vmin <<"\t" << size <<std::endl;
  
    for (GLuint i=0;i<size;i++){
      color[3*i]=color[3*i+1]=color[3*i+2]=1.0;
      double v=s[i];
     
      if (v < (vmin + 0.25 * dv)) {
	color[3*i] = 0;
	color[3*i+1] = 4 * (v - vmin) / dv;
      } else if (v < (vmin + 0.5 * dv)) {
	color[3*i] = 0;
	color[3*i+2] = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
      } else if (v < (vmin + 0.75 * dv)) {
	color[3*i] = 4 * (v - vmin - 0.5 * dv) / dv;
	color[3*i+2] = 0;
      } else {
	color[3*i+1] = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
	color[3*i+2] = 0;
      }
     
    }
  
  }

  MeshRamper::MeshRamper(TriMesh *pmesh, ScalarFunction* psfun)
    :MeshPainter(pmesh){
    color_array = new GLfloat[3*vn];
    color_ramping(color_array, psfun);  
  };

  /*
  MeshRamper::MeshRamper(TriMesh *pmesh, ScalarFunction* psfun, bool equalhist)
    :MeshPainter(pmesh){
    color_array = new GLfloat[3*vn];
    color_ramping(color_array, psfun, equalhist);  
  };
  */


  void MeshRamper::prepare(){
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    //load geometric information from painte

    glNormalPointer(GL_DOUBLE, 0, normal_array);
    glVertexPointer(3, GL_DOUBLE, 0, vertex_array); 
    glColorPointer(3, GL_FLOAT, 0, color_array);

  }



  MeshRamper::~MeshRamper() {
    delete [] color_array;
  }

  MeshMarker::MeshMarker(TriMesh *pmesh, BooleanFunction* pbfun)
    :MeshPainter(pmesh) {
    mark_count = std::accumulate(pbfun->begin(), pbfun->end(), 0);
    mark_array = new GLuint[mark_count];    
    int j=0, n=pbfun->size();
    for (int i=0; i < n; ++i) 
      if ((*pbfun)[i]) mark_array[j++]=i;
    
  }

  MeshMarker::MeshMarker(TriMesh *pmesh, std::vector<int> &index_mark)
    :MeshPainter(pmesh) {
    mark_count = index_mark.size();
    mark_array = new GLuint[mark_count];    
    for (unsigned i=0; i<mark_count; ++i) mark_array[i] = index_mark[i];
  }

  void MeshMarker::prepare(){
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

    glNormalPointer(GL_DOUBLE, 0, normal_array);
    glVertexPointer(3, GL_DOUBLE, 0, vertex_array); 
  }

  void MeshMarker::draw(){

    glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.3, 0.5, 0.6);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glCullFace(GL_BACK);
  

    glDrawElements(GL_TRIANGLES, 3*fn, GL_UNSIGNED_INT, index_array);

    //glDisable(GL_LIGHTING);
    glColor3f(1., 0., 0.);
    //glPointSize(4);
    //glDrawElements(GL_POINTS, mark_count, GL_UNSIGNED_INT, mark_array);

    GLfloat radius = 2 * (MeshViewer::currentMeshViewer->coordinate_max_x - MeshViewer::currentMeshViewer->coordinate_min_x)/ MeshViewer::currentMeshViewer->width / scale * perfect_factor;
    
    for (GLuint i= 0; i< mark_count; ++i){
      glPushMatrix();
      glTranslatef(vertex_array[3*mark_array[i]], vertex_array[3*mark_array[i] +1] , vertex_array[3*mark_array[i]+2]);
      glutSolidSphere(radius , 20, 20);
      glPopMatrix();
    }
    glDisable(GL_COLOR_MATERIAL);

    //glEnable(GL_LIGHTING);

    
    
  }
  
  MeshMarker::~MeshMarker() {
    delete [] mark_array;
  }





  MeshViewer::MeshViewer(int argc, char** argv) 
    :width(800), height(800) {
    currentMeshViewer = this; //static member need definition

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowPosition(200.0, 0.0);
    glutInitWindowSize(width, height);
    glutCreateWindow("tri-mesh-toolkit - Viewer");

    glShadeModel(GL_SMOOTH);// Enable Smooth Shading
    
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
    GLuint n = currentMeshViewer->Painters.size();
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers

    glPushMatrix();
    
    for (GLuint i=0; i<n; i++) 
      { 

	glPushMatrix();//push i-th matrix
	glCallList(currentMeshViewer->Painters[i]->LIST_NAME);

      
	//currentMeshViewer->Painters[i]->prepare();    
	//currentMeshViewer->Painters[i]->draw();
	glPushMatrix();//pop i-th matrix

      }

    glPopMatrix();

    glutSwapBuffers();

    
      
  }

  
  void MeshViewer::reshape(int width, int height){

    glViewport(0, 0, width, height);    
    
    GLfloat radio_x = (currentMeshViewer->coordinate_max_x - currentMeshViewer->coordinate_min_x) /  width;
    GLfloat radio_y = (currentMeshViewer->coordinate_max_y - currentMeshViewer->coordinate_min_y) /  height;
    GLfloat radio = scale *((radio_x > radio_y)? radio_x: radio_y);




    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho( currentMeshViewer->center_x - radio * width/perfect_factor, 
	     currentMeshViewer->center_x + radio * width/perfect_factor,
	     currentMeshViewer->center_y - radio * height/perfect_factor, 
	     currentMeshViewer->center_y + radio * height/perfect_factor,
	     //100000, -100000);
	     - currentMeshViewer->center_z - 2* currentMeshViewer->length_z ,
	     - currentMeshViewer->center_z + 2* currentMeshViewer->length_z);//(NEW) set up our viewing area


    //glOrtho(-1,1,-1,1,-10,10);
    glMatrixMode(GL_MODELVIEW);

    
  }

  void MeshViewer::keyboard(unsigned char key, int x, int y) {
    if (key=='q'||key=='Q') exit(0);
  }



  void MeshViewer::motion(int x, int y) {
    if (left_button) {
      glLoadIdentity();
      glTranslated(
		   (currentMeshViewer->coordinate_max_x - currentMeshViewer->coordinate_min_x) *(GLfloat) (x-origin_x)/(GLfloat) currentMeshViewer->width * scale * 2./perfect_factor,
		    - (currentMeshViewer->coordinate_max_y - currentMeshViewer->coordinate_min_y) * (GLfloat) (y-origin_y)/(GLfloat) currentMeshViewer->height * scale* 2./ perfect_factor, 0.0);
      glMultMatrixf(CTM);

      glutPostRedisplay();
    }

    else if (right_button) {
      glLoadIdentity();

      glRotatef(90.0* (GLfloat) (x-origin_x)/(GLfloat) currentMeshViewer->width * scale * 2./perfect_factor, 0.0, 1.0, 0.0);
      glRotatef(90.0* (GLfloat) (y-origin_y)/(GLfloat) currentMeshViewer->height *scale * 2./perfect_factor, 1.0, 0.0, 0.0);                

      glMultMatrixf(CTM);	

      glutPostRedisplay();
    }

  }


// compatibility with original GLUT

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

#define scale_coeff       (1.2)

  void MeshViewer::mouse(int button, int state, int x, int y) {
    switch(button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {
	origin_x = x; origin_y = y; left_button =true;
	glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      }
      else if (state == GLUT_UP) {
	left_button = false;
	glutPostRedisplay();
      }
      break;
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
	glLoadIdentity();
	scale =1.;
	reshape(currentMeshViewer->width, currentMeshViewer->height);
	glutPostRedisplay();
      }
      break;
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
	origin_x = x; origin_y = y; right_button =true;
	glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      }
      else if (state == GLUT_UP) {
	right_button = false;
	glutPostRedisplay();
      }

      break;
    case GLUT_WHEEL_UP:
      if (state == GLUT_UP){	
	scale *= scale_coeff;
	reshape(currentMeshViewer->width, currentMeshViewer->height);
	glutPostRedisplay();      
      }
      break;
    case GLUT_WHEEL_DOWN:
      if (state == GLUT_DOWN){	
	scale /= scale_coeff;
	reshape(currentMeshViewer->width, currentMeshViewer->height);
	glutPostRedisplay();      
      }
      break;
    default:
      break;
    }

  }



  void MeshViewer::init(){

  
    /* 3D configuration */
    // clear the window color and depth buffer
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClearDepth(1.0);
    
    center_x = (coordinate_min_x + coordinate_max_x) /2.;
    center_y = (coordinate_min_y + coordinate_max_y) /2.;
    center_z = (coordinate_min_z + coordinate_max_z) /2.;
    length_z = (coordinate_max_z - coordinate_min_z) /2.;
  
    GLfloat radio_x = (coordinate_max_x - coordinate_min_x) /  width;
    GLfloat radio_y = (coordinate_max_y - coordinate_min_y) /  height;
    GLfloat radio = (radio_x > radio_y)? radio_x: radio_y;

    //std::cout<< length_z << std::endl;


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho( center_x - radio * width/perfect_factor, center_x + radio * width/perfect_factor,
	     center_y - radio * height/perfect_factor, center_y + radio * height/perfect_factor,
	     //100000, -100000);
	     - center_z - 2* length_z , - center_z + 2* length_z);//(NEW) set up our viewing area


    //glOrtho(-1,1,-1,1,-10,10);
    glMatrixMode(GL_MODELVIEW);


    /////////////////////////
    scale =1.;


    ///////////////////////
    add_lights();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);

  }


  void MeshViewer::add_lights(){
    /*
      light_position0[0] = coordinate_max_x * perfect_factor;
      light_position0[1] = coordinate_max_y * perfect_factor;
      light_position0[2] = coordinate_max_z * perfect_factor;

      light_position1[0] =  coordinate_min_x * perfect_factor;
      light_position1[1] =  coordinate_min_y * perfect_factor;
      light_position1[2] =  coordinate_min_z * perfect_factor;
    */

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);


  }
  void MeshViewer::add_painter(MeshPainter *painter){
    Painters.push_back(painter);
    GLuint count = Painters.size();
    if (painter->coordinate_min_x < coordinate_min_x|| count ==1) 
      coordinate_min_x = painter->coordinate_min_x;
    if (painter->coordinate_min_y < coordinate_min_y|| count ==1) 
      coordinate_min_y = painter->coordinate_min_y;
    if (painter->coordinate_min_z < coordinate_min_z|| count ==1)
      coordinate_min_z = painter->coordinate_min_z;
    if (painter->coordinate_max_x > coordinate_max_x|| count ==1) 
      coordinate_max_x = painter->coordinate_max_x;
    if (painter->coordinate_max_y > coordinate_max_y|| count ==1) 
      coordinate_max_y = painter->coordinate_max_y;
    if (painter->coordinate_max_z > coordinate_max_z|| count ==1) 
      coordinate_max_z = painter->coordinate_max_z;


    painter->prepare();    

    glNewList(painter->LIST_NAME, GL_COMPILE_AND_EXECUTE);
    painter->draw();
    glEndList();


  }

  void MeshViewer::view(){

    glutMainLoop();
  }













}
