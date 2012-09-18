/*
  FILE: ManifoldTriMesh.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implemented class ManifoldTriMesh.
  
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

#include "meshtk/ManifoldTriMesh.hh"
namespace meshtk {
  void ManifoldTriMesh::init_index(){
    TriMesh::init_index();
    
    vertex_coord_seq = PointFuncSequence(1);
    vertex_coord_seq[0] = PointFunction(vertex_num);

    for (int i=0;i<vertex_num; i++){
      vertex_coord_seq[0][i] =  IV[i]->point();
    }

  }
  

}
