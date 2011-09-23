/*
  FILE: DynamicTriMesh.hh This file is part of MeshTK.
  It is a C++ source file which implemented class DynamicTriMesh, achieving smoothing, 
  denoising and remeshing.
  
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

#include "meshtk/DynamicTriMesh.hh"



namespace meshtk {
  
  void DynamicTriMesh::init_index(){
    TriMesh::init_index();
    
    vertex_coord = PointFunction(vertex_num);
    facet_coord = PointFunction(facet_num);
  }

  
  void DynamicTriMesh::restore_coord(){
    for (int i=0;i<vertex_num; i++){
      IV[i]->point() = vertex_coord[i];
    }
  }

  void DynamicTriMesh::gaussian_smooth(double coeff){
  }

}
