/*
  FILE: ex0.cc This file is part of tri-mesh-toolkit.
  This example shows how to manipulate some basic functions of libmeshtk

  ./ex0 mesh out

  will read mesh.off and write to out.off
  
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


// include the main header file 
#include <meshtk/TriMesh.hh>


int main(int argc, char *argv[])
{
  // define a mesh object with empty content
  meshtk::TriMesh Mesh;

  // read mesh from file specified by name and type
  Mesh.read(argv[1], "off");// load mesh.off into Mesh 

  // compulsory initialization, with initializing the inner index system
  // (for vertex, facet and halfedge). But it does not update geometric information
  Mesh.init_index(); // should be called before any other mesh processing routines
  
  // update some basic geometric information
  Mesh.update_base();
  // output some basic information
  std::cout << "number of vertices: "<< Mesh.vertex_num << std::endl;
  std::cout << "number of facets: "<< Mesh.facet_num << std::endl;
  std::cout << "number of halfedges: "<< Mesh.halfedge_num << std::endl;
  std::cout << "total area: "<< Mesh.total_area << std::endl;

  // print vertex normal computed by update_base() to file
  Mesh.attribute_print(MESHTK_VERTEX_NORM, MESHTK_VECTOR, argv[2]);

  // write mesh to file
  Mesh.write(argv[2], "off");
  
  return 0;
}




