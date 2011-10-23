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
  double TriMesh::update_vertex_geodesic_distance(int source_vertex_index,
						  ScalarFunction & geodesic_distance){

    clock_start("Compute vertex-source geodesic distance");

    int vertex_size = geodesic_distance.size();
    // pre condition with function update_compact_base()
    geodesic::Mesh mesh; 
    mesh.initialize_mesh_data(vertex_array, tri_index_array);
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);

    geodesic::SurfacePoint source(&mesh.vertices()[source_vertex_index]); //create source 
    std::vector<geodesic::SurfacePoint> single_source(1,source); 

    // propagation with certain distance
    algorithm.propagate(single_source);	

    double avg_geodesic_distance = 0;
    for (int i=0; i < vertex_size; ++i) {
      geodesic::SurfacePoint p(&mesh.vertices()[i]);

      //for a given surface point, find closets source and distance to this source
      algorithm.best_source(p,geodesic_distance[i]);		
      avg_geodesic_distance += geodesic_distance[i];
    }

    clock_end();
    return avg_geodesic_distance / vertex_size;
  }

    
  int TriMesh::assemble_export_Nystrom_matrix(std::vector<int> & sampling,
					      int addition_size,
					      std::string name,
					      double (*distance_function) (int, ScalarFunction&)) {
    int init_size = sampling.size();
    int total_size = init_size;

    double max_B_distance = 0;
    int max_B_distance_index =  vertex_num * rand()/ RAND_MAX;

    ScalarFunction distance;
    ScalarFunction new_distance;

    distance.resize(vertex_num);
    for (int i =0;i<vertex_num; ++i) distance[i] =-1;

    std::string matdata_name = name; 
    std::ofstream mat_data(matdata_name.c_str(), std::ios::out|std::ios::binary);

    std::string matinfo_name = name; matinfo_name.append(".info");
    std::ofstream mat_info(matinfo_name.c_str());

    mat_info << "#size " << vertex_num << "\n";

    for (int i = 0; i < init_size; ++i ){
      new_distance.clear(); 
      new_distance.resize(vertex_num); 

      double distance_avg = (*distance_function) (sampling[i], new_distance);


      mat_info << sampling[i] <<"\n";

      for (int j=0; j< vertex_num; ++j) 
	if (new_distance[j] < distance[j] || distance[j] <0) {
	  distance[j] = new_distance[j];
	}

      for (int j=0; j< vertex_num; ++j) {
	new_distance[j] *= std::sqrt( facet_area[j] * facet_area[sampling[i]]);

      }
      mat_data.write((char *) &new_distance[0], vertex_num * sizeof(MESHTK_SCALAR_TYPE));
    }

    for (int j=0; j<vertex_num; ++j) 
      if (distance[j] > max_B_distance) {
	max_B_distance = distance[j];
	max_B_distance_index = j;
      }

    

    while (total_size < (init_size + addition_size)){
      new_distance.clear(); //new_SIFT_distance.clear();
      new_distance.resize(vertex_num); 
      //new_SIFT_distance.resize(vertex_num);

      sampling.push_back(max_B_distance_index);

      double distance_avg = (*distance_function) (max_B_distance_index, new_distance);

      mat_info << max_B_distance_index <<"\n";

      for (int j=0; j< vertex_num; ++j) 
	if (new_distance[j] < distance[j] || distance[j] <0) {
	  distance[j] = new_distance[j];
	}
      
      for (int j=0; j< vertex_num; ++j) {
	new_distance[j] *= std::sqrt( facet_area[j] * facet_area[max_B_distance_index]);
      }
      mat_data.write((char *) &new_distance[0], vertex_num * sizeof(MESHTK_SCALAR_TYPE));
      
      max_B_distance = 0;
      for (int j=0; j<vertex_num; ++j) 
	if (distance[j] > max_B_distance) {
	  max_B_distance = distance[j];
	  max_B_distance_index = j;
	}
      ++ total_size;

    }

    mat_info.close();
    mat_data.close();


    return total_size;
  }


}



