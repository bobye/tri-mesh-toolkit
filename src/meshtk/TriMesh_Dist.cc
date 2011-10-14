/*
  FILE: TriMesh_Curv.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implements the distances between vertices on mesh
  of class TriMesh.
  
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


#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"
// include header file for geodesic computation
#include <geodesic/geodesic_algorithm_exact.h>


namespace meshtk {
  double TriMesh::update_vertex_geodesic(int source_vertex_index,
				ScalarFunction & geodesic_distance){

    clock_start("Compute vertex-source geodesic distance");


    // pre condition with function update_compact_base()
    geodesic::Mesh mesh; 
    mesh.initialize_mesh_data(vertex_array, tri_index_array);
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);

    geodesic::SurfacePoint source(&mesh.vertices()[source_vertex_index]); //create source 
    std::vector<geodesic::SurfacePoint> single_source(1,source); 

    // propagation with certain distance
    algorithm.propagate(single_source);	

    double avg_geodesic_distance = 0;
    for (int i=0; i < vertex_num; ++i) {
      geodesic::SurfacePoint p(&mesh.vertices()[i]);

      //for a given surface point, find closets source and distance to this source
      algorithm.best_source(p,geodesic_distance[i]);		
      avg_geodesic_distance += vertex_area[i] * geodesic_distance[i];
    }

    clock_end();
    return avg_geodesic_distance / total_area;
  }

  double TriMesh::update_vertex_biharmonic(int source_vertex_index,
					   ScalarFunction & biharmonic_distance) {
    
    return 0;
  }
}



