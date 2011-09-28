/*
  FILE: test.cc This file is part of MeshTK.
  It is a C++ source file which implements some experiemtal functions.
  functions of meshtk library.
  
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

#include <string>
#include <algorithm>

#include <meshtk/TriMesh.hh>
#include <meshtk/DynamicTriMesh.hh>
#include <meshtk/MeshViewer.hh>

int main(int argc, char *argv[])
{
  /***************************************************************************/    
  // region to test
  //mesh.update_vertex_salient(5,1);
    
  meshtk::DynamicTriMesh mesh;
  meshtk::DynamicTriMesh mesh_base;
  mesh.read("mesh", "off");
  mesh_base.read("mesh_base", "off");

  mesh.init_index();
  mesh_base.init_index();
  mesh.update_base();
  mesh_base.update_base();
  mesh.load_coord();
  mesh_base.load_coord();

  meshtk::PointFunction *mesh_coord = (meshtk::PointFunction *) mesh.attribute_extract(MESHTK_VERTEX_COORD);
  meshtk::PointFunction *mesh_base_coord=(meshtk::PointFunction *) mesh_base.attribute_extract(MESHTK_VERTEX_COORD);

  unsigned USER_MESH_OFFSET = mesh_base.attribute_allocate(MESHTK_VERTEX, MESHTK_VECTOR);
  meshtk::VectorFunction *mesh_offset = (meshtk::VectorFunction *) mesh_base.attribute_extract(USER_MESH_OFFSET);

  unsigned USER_MESH_OFFSET_NORM = mesh_base.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
  meshtk::ScalarFunction *mesh_offset_proj = (meshtk::ScalarFunction *) mesh_base.attribute_extract(USER_MESH_OFFSET_NORM);

  meshtk::VectorFunction *mesh_base_norm = (meshtk::VectorFunction *) mesh_base.attribute_extract(MESHTK_VERTEX_NORM);

  int vertex_size = mesh_coord->size();

  for (int i=0;i<vertex_size;++i) {
    (*mesh_offset)[i] = (*mesh_coord)[i] - (*mesh_base_coord)[i];
    (*mesh_offset_proj)[i] = (*mesh_offset)[i] * (*mesh_base_norm)[i];
    //(*mesh_base_coord)[i] = (*mesh_base_coord)[i] + (*mesh_offset_proj)[i] * (*mesh_base_norm)[i];
  }

  unsigned USER_MESH_BUFFER_SCALAR = mesh_base.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
  meshtk::ScalarFunction *mesh_buffer_scalar = (meshtk::ScalarFunction *) mesh_base.attribute_extract(USER_MESH_BUFFER_SCALAR);

  // smooth iteration
  mesh_base.update_vertex_neighbor(3.);
  for (int i=0; i<10; i++) {
    std::cout<<i<<std::endl;
    mesh_base.gaussian_smooth_vertex(1., *mesh_offset_proj, *mesh_buffer_scalar, 0.);
    std::swap(mesh_offset_proj, mesh_buffer_scalar);
  }
  //mesh_base.restore_coord();
  //mesh.write("out","off");


  /*  
  for (int i=0;i<vertex_size;++i) (*mesh_base_coord)[i] = Point(0,0,0)+ ((*mesh_coord)[i] - (*mesh_base_coord)[i]);
  mesh_base.restore_coord();
  mesh_base.update_base();

  mesh_base.write("out","off");
  */

  
  meshtk::MeshViewer viewer(argc, argv);
  meshtk::MeshRamper painter(&mesh_base, mesh_offset_proj);
  //meshtk::MeshPainter painter(&mesh_base);
  viewer.add_painter(&painter);

  viewer.init();// call this func last before loop
  viewer.view();
    
    


  return 0;
}
