/*
  FILE: TriMesh_Curv.cc This file is part of MeshTK.
  It is a C++ source file which implements the local point feature 
  of class TriMesh.
  
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

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"

// include header file for geodesic computation
#include <geodesic/geodesic_algorithm_exact.h>

namespace meshtk{

  struct front_circ {
    Halfedge_handle hf; // Halfedge 
    front_circ *prev; // pointers
    front_circ *next;
    bool fr;// front reached?
  };

  int TriMesh::update_vertex_neighbor_euclidean(int source_vertex_index, 
						double propagation_distance,
						ScalarNeighborFunction &vertex_neighbor,
						NeighborIndex &facet_neighbor){
    
    //NeighborIndex vertex_outbound;
    Vector displace; double distance;
    front_circ *start = new front_circ; //start->prev = start->next = start;
    front_circ *front = start;
    vertex_neighbor[source_vertex_index] = 0.;

    HV_circulator hv = IV[source_vertex_index]->vertex_begin();
    do {	
      front->next = new front_circ;
      front->next->prev = front; 
      front = front->next;
      front->hf = hv->prev()->opposite();

      //if (hv->facet()!=NULL) facet_neighbor_euclidean[i].insert(hv->facet()->index);
      int vertex_index = hv->opposite()->vertex()->index;
      displace = IV[vertex_index]->point() - IV[source_vertex_index]->point();
      distance = std::sqrt(displace * displace);

      vertex_neighbor[vertex_index] = distance;

      if (hv->facet()!=NULL ) facet_neighbor.insert(hv->facet()->index);

    } while (++hv != IV[source_vertex_index]->vertex_begin());
    start->next->prev = front; front->next = start->next;
    delete start;
    start = front;
    bool front_reached = false;
      
    do {	  
      //if (i==M) getchar();
	
      front_reached = true;


      while (front->hf->facet()==NULL){
	// delete verbose boundary halfedge to improve performance and 
	// insert the corresponding vertex into neighbor set
	front = front->next;
	front->prev = front->prev->prev;
	delete front->prev->next;
	front->prev->next = front;
	start = front;
      }

      // to test whether propagate toward the facet
      int facet_index, vertex_index = front->hf->next()->vertex()->index;
      if (facet_neighbor.find(facet_index = front->hf->facet()->index) == facet_neighbor.end() && !front->fr) {
	// && vertex_outbound.find(vertex_index = front->hf->next()->vertex()->index) == vertex_outbound.end()) {

	displace = IV[vertex_index]->point() - IV[source_vertex_index]->point();
	distance = std::sqrt(displace * displace);
	if (distance < propagation_distance) {
	  // insert new vertex to the front
	  facet_neighbor.insert(facet_index);
	  vertex_neighbor[vertex_index] = distance;
	    
	  front_circ *insert = new front_circ;
	  insert->hf = front->hf->prev()->opposite();
	  insert->prev = front->prev; insert->next = front;
	  front->hf = front->hf->next()->opposite();
	  front->prev->next = insert; front->prev = insert;
	  front_reached =false;

	  start=front;
	  //if (i==M) std::cout<<"+ " <<front->hf->opposite()->vertex()->index << " -> " <<front->hf->vertex()->index << std::endl;
	}
	else {
	  front->fr = true; 
	  front = front->next;
	}

      }	

      else {
	front = front->next; 
      }

      while (front->hf->opposite() == front->next->hf) {	  
	// delete verbose pair of halfedges to improve performance and 
	// insert the corresponding vertex into neighbor set
	//if (i==M) std::cout<<"- "<<front->hf->vertex()->index <<std::endl;
	  
	front->prev->next = front->next->next;
	front = front->prev;
	delete front->next->prev->prev;
	delete front->next->prev;
	front->next->prev = front;
	start = front;
      }

    } while (front!= start || !front_reached);
      
    front->next->prev = NULL;
    while (front->prev!= NULL) { front = front->prev; delete front->next;  } delete front;
    return vertex_neighbor.size();

  }


  double TriMesh::update_vertex_neighbor_euclidean(double propagation_distance_coeff){

    clock_start("Update neighbors w.r.t Euclidean");


    double distance_threshold = propagation_distance_coeff * avg_edge_len; 

    //vertex_neighbor.clear();
    //vertex_neighbor.resize(vertex_num);
    facet_neighbor_euclidean.clear();
    facet_neighbor_euclidean.resize(vertex_num);
    vertex_neighbor_euclidean.clear();
    vertex_neighbor_euclidean.resize(vertex_num);
    int count=0;// to record average neighbor size of vertices

    // The algorithm implemented below is a wide search 
    // which is not very efficient
    /*
    for (int i = 0; i < vertex_num; ++i) {
      NeighborIndex buffer;      
      buffer.insert(i);
      while (!buffer.empty()) {
	NeighborIndex::iterator it = buffer.begin();
	Vector displace = IV[*it]->point()- IV[i]->point();
	
	if (displace * displace < distance_threshold * distance_threshold) {
	  vertex_neighbor[i].insert(*it);
	  
	  HV_circulator hv = IV[*it]->vertex_begin();
	  do {
	    int j =hv->opposite()->vertex()->index;
	    if (vertex_neighbor[i].find(j) == vertex_neighbor[i].end()) buffer.insert(j);
	  }while (++hv != IV[*it]->vertex_begin());
	}

	buffer.erase(it);	
      }

      count += vertex_neighbor[i].size();
    }
    */

    // the algorithm implemented below is a front propagation one
    // initialization from a 1-ring neighbor

    

    for (int i = 0; i < vertex_num; ++i){      
      count += update_vertex_neighbor_euclidean(i, 
						distance_threshold,
						vertex_neighbor_euclidean[i],
						facet_neighbor_euclidean[i]);
    }


    clock_end();
    return (double) count/(double) vertex_num;
  }


  /*
  double TriMesh::local_quadratic_extrema(ScalarFunction &value, int focus){
    return 0;
  }
  */

  /*
  void TriMesh::geodesic_init() {
    //create internal mesh data structure including edges
    geodesic_mesh = new geodesic::Mesh;
    geodesic_mesh -> initialize_mesh_data(vertex_array, tri_index_array);		
    //create exact algorithm for the mesh
    geodesic_algorithm = new geodesic::GeodesicAlgorithmExact(geodesic_mesh);	  
  }
  */


  int TriMesh::update_vertex_neighbor_geodesic(int source_vertex_index, 
					       double propagation_distance,
					       ScalarNeighborFunction &vertex_neighbor,
					       NeighborIndex &facet_neighbor,
					       ScalarNeighborFunction &vertex_neighbor_interest,
					       NeighborIndex &facet_neighbor_interest) {

    std::vector<double> points(3*vertex_neighbor_interest.size());
    std::vector<unsigned> faces(3*facet_neighbor_interest.size());
    std::map<unsigned, unsigned> local_index;
    int local_source_vertex_index, j=0;

    for (ScalarNeighborFunction::iterator it=vertex_neighbor_interest.begin();
	 it != vertex_neighbor_interest.end(); ++it, ++j) {

      points[3*j] = vertex_array[3 * (it->first)];
      points[3*j+1] = vertex_array[3 * (it->first)+1];
      points[3*j+2] = vertex_array[3 * (it->first)+2];
      local_index[it->first] = j;
      if (it->first == source_vertex_index) local_source_vertex_index = j;

    }

    j=0;
    for (NeighborIndex::iterator it=facet_neighbor_interest.begin();
	 it!=facet_neighbor_interest.end(); ++it, ++j) {

      faces[3*j] = local_index[tri_index_array[3* (*it)]];
      faces[3*j+1] = local_index[tri_index_array[3* (*it) +1]];
      faces[3*j+2] = local_index[tri_index_array[3* (*it) +2]];

    }

  
    geodesic::Mesh mesh; 
    mesh.initialize_mesh_data(points, faces);
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);

    geodesic::SurfacePoint source(&mesh.vertices()[local_source_vertex_index]); //create source 
    //in general, there could be multiple sources, but now we have only one
    std::vector<geodesic::SurfacePoint> single_source(1,source); 

    // propagation with certain distance
    algorithm.propagate(single_source, propagation_distance);	

    facet_neighbor = facet_neighbor_interest; // copy the elements.

    for(ScalarNeighborFunction::iterator it = vertex_neighbor_interest.begin();
	it != vertex_neighbor_interest.end(); ++it) {
      geodesic::SurfacePoint p(&mesh.vertices()[local_index[it->first]]);		
      
      double distance;
      algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source
      if (distance < propagation_distance) // it is within the geodesic range
	vertex_neighbor[it->first] = distance;
      else { // otherwise remove its associated facets from facet_neighbor
	HV_circulator hv = IV[it->first]->vertex_begin();
	do {
	  if (hv->facet()!=NULL) facet_neighbor.erase(hv->facet()->index);
	} while (++hv != IV[it->first]->vertex_begin());
      }
    }

    return vertex_neighbor.size();

  }


  double TriMesh::update_vertex_neighbor_geodesic(double propagation_distance_coeff) {
    // it is very slow to compute a neighbor with geodesic threshold

    clock_start("Update neighbors w.r.t geodesic"); 

    double propagation_distance = propagation_distance_coeff * avg_edge_len;
    int count = 0;
    vertex_neighbor_geodesic.clear();
    vertex_neighbor_geodesic.resize(vertex_num);


    for (int i=0;i < vertex_num; ++i) {
      //printf("\b\b\b%2d%%", (i+1)*100 / vertex_num);

      count += update_vertex_neighbor_geodesic(i, 
					       propagation_distance, 
					       vertex_neighbor_geodesic[i],
					       facet_neighbor_geodesic[i],
					       vertex_neighbor_euclidean[i],
					       facet_neighbor_euclidean[i]);
      //std::cout << i << std::endl;
    }

    clock_end();

    return (double) count / (double) vertex_num;
  }

  void TriMesh::register_vertex_keypoint(int vertex_index, 
				double scale_distance,
				double magnitude,
				ScalarFunction& scale_space_function) {
    ScalarNeighborFunction buffer_vertex_neighbor_euclidean;
    NeighborIndex buffer_facet_neighbor_euclidean;
    
    keypoints.push_back(KeyPoint(vertex_index, scale_distance, magnitude));

    update_vertex_neighbor_euclidean(vertex_index, 
				     scale_distance,
				     buffer_vertex_neighbor_euclidean,
				     buffer_facet_neighbor_euclidean);

    update_vertex_neighbor_geodesic(vertex_index,
				    scale_distance,
				    keypoints.back().vertex_neighbor,
				    keypoints.back().facet_neighbor,
				    buffer_vertex_neighbor_euclidean,
				    buffer_facet_neighbor_euclidean);
    
  }

  int TriMesh::detect_vertex_keypoint(ScalarFunction &valueScalar, 
				      BooleanFunction &keyBoolean, 
				      int iter, 
				      int pre_iter){
    // in the implementation of keypoint detection, no threshold is taken.
    int count = 0;
    double coeff = .6;
    double neighbor_size = 3*coeff > 3.? 3*coeff: 3.;

    update_vertex_neighbor_euclidean(neighbor_size);     
    update_vertex_neighbor_geodesic(neighbor_size); 
    
    neighbor_distance_map  = & vertex_neighbor_geodesic;
    //neighbor_distance_map  = & vertex_neighbor_euclidean;


    clock_start("Start keypoint detection");


    // update static coefficient for smooth
    double sigma = coeff * avg_edge_len;
    std::vector<ScalarNeighborFunction > coeff_list(vertex_num);


    for (int i=0; i<vertex_num; ++i) {
      double scale, total_scale = 0;

      
      for (ScalarNeighborFunction::iterator it = (*neighbor_distance_map)[i].begin();
	   it != (*neighbor_distance_map)[i].end(); ++it) {
	// use the euclidean distance to contruct smooth stencil
	scale = it->second;
	scale = vertex_area[it->first] * std::exp( - (scale * scale) / (2 * sigma * sigma));
	total_scale += scale;
	coeff_list[i][it->first] = scale;
      }

      for (ScalarNeighborFunction::iterator it = coeff_list[i].begin();
	   it != coeff_list[i].end(); ++it) {
	it->second /= total_scale;
      }
	
    }



    ScalarFunction buffer_v[4];
    ScalarFunction buffer_dv[4];

    for (int i=0; i<4; ++i) {buffer_v[i].resize(vertex_num); buffer_dv[i].resize(vertex_num);}
    int buffer_curr = 0;
    
    buffer_v[0] = valueScalar;

    /*
    gaussian_smooth_vertex(coeff, buffer_v[0], buffer_v[1], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[1], buffer_v[2], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[2], buffer_v[3], 0.);
    */
    
    smooth_vertex(coeff_list, buffer_v[0], buffer_v[1], 0.);
    smooth_vertex(coeff_list, buffer_v[1], buffer_v[2], 0.);
    smooth_vertex(coeff_list, buffer_v[2], buffer_v[3], 0.);


    for (int i = 0; i < vertex_num; ++i) {
      buffer_dv[0][i] =   (buffer_v[1][i] - buffer_v[0][i]) / std::log(2.);
      buffer_dv[1][i] =   (buffer_v[2][i] - buffer_v[1][i]) / std::log(1.5);
    }
    
    for (int i=0; i< iter; ++i) {
      
      //      double radio = 0,// (0.001 / std::sqrt(total_area)) * coeff * coeff, 
      double radio2 = 0.05 * coeff * coeff;

      int zero = buffer_curr, one=(buffer_curr +1)%4, 
	two=(buffer_curr +2)%4, three=(buffer_curr +3)%4;

      //gaussian_smooth_vertex(coeff, buffer_v[two], buffer_v[three], 0.);
      smooth_vertex(coeff_list, buffer_v[two], buffer_v[three], 0.);


      for (int j=0; j < vertex_num; ++j) 
	buffer_dv[two][j] = (buffer_v[three][j] - buffer_v[two][j]) / std::log(double (i+3)/ double(i+2));




      if (pre_iter--<=0 && i*coeff >= 1.) {
	for (int j=0; j < vertex_num; ++j) {
	    HV_circulator hv=IV[j]->vertex_begin();
	    bool tmp = true;
	  
	    if (buffer_dv[one][j] < buffer_dv[two][j]  || buffer_dv[one][j] < buffer_dv[zero][j] )
	      tmp = false;	  
	    else do {
		if (buffer_dv[one][j] < buffer_dv[one][hv->next()->vertex()->index] ||
		    buffer_dv[one][j] < buffer_dv[two][hv->next()->vertex()->index] ||
		    buffer_dv[one][j] < buffer_dv[zero][hv->next()->vertex()->index]) {
		  tmp = false; break;
		}
	      }while (++hv!=IV[j]->vertex_begin());
	    if (tmp && std::fabs(buffer_dv[one][j]) * avg_edge_len > radio2 )
	      { keyBoolean[j] = true; ++ count; 
		register_vertex_keypoint(j, 
					 std::sqrt(i+2)* sigma,
					 std::fabs(buffer_dv[one][j] * avg_edge_len / radio2),
					 buffer_v[two]);
		continue;}

	  
	    hv=IV[j]->vertex_begin(); tmp = true;
	    if (buffer_dv[one][j] > buffer_dv[two][j]  || buffer_dv[one][j] > buffer_dv[zero][j] )
	      tmp = false;
	    else do {
		if (buffer_dv[one][j] > buffer_dv[one][hv->next()->vertex()->index] ||
		    buffer_dv[one][j] > buffer_dv[two][hv->next()->vertex()->index] ||
		    buffer_dv[one][j] > buffer_dv[zero][hv->next()->vertex()->index]) {
		  tmp = false; break;
		}
	      }while (++hv!=IV[j]->vertex_begin());
	    if (tmp && std::fabs(buffer_dv[one][j]) * avg_edge_len  > radio2 )
	      { keyBoolean[j] = true; ++ count; 
		register_vertex_keypoint(j, 
					 std::sqrt(i+2)* sigma,
					 std::fabs(buffer_dv[one][j]) * avg_edge_len / radio2,
					 buffer_v[two]);
		continue; }

	  }
	
	std::cout<< "#" << i << " detect Keypoint iteration\t Found " << count << std::endl;
	
      } else {
	std::cout <<"#" << i << " pre-smooth count" << std::endl;
      }      

      buffer_curr = (buffer_curr +1)%4;

    }

    valueScalar = buffer_v[ buffer_curr ];

    clock_end();
    return count;
  }

}


