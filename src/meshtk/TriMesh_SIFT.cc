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
#include <ctime>

namespace meshtk{


  double TriMesh::update_vertex_neighbor(double coeff){

    time_t  start, end; time(&start);
    std::cout << "Update neighbors ... " << std::flush; 

    double square_distance_threshold = coeff * avg_edge_len; 
    square_distance_threshold *=square_distance_threshold;

    vertex_neighbor.clear();
    vertex_neighbor.resize(vertex_num);
    int count=0;// to record average neighbor size of vertices

    // The algorithm implemented below is a wide search 
    // which is not very efficient
    /*
    for (int i = 0; i < vertex_num; ++i) {
      std::set<int> buffer;      
      buffer.insert(i);
      while (!buffer.empty()) {
	std::set<int>::iterator it = buffer.begin();
	Vector displace = IV[*it]->point()- IV[i]->point();
	
	if (displace * displace < square_distance_threshold) {
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

    struct front_circ {
      Halfedge_handle hf; // Halfedge 
      front_circ *prev; // pointers
      front_circ *next;
      bool fr;// front reached?
    };
    

    for (int i = 0; i < vertex_num; ++i){      
      std::set<int> facet_occupied;
      //std::set<int> vertex_outbound;
      Vector displace;
      front_circ *start = new front_circ; //start->prev = start->next = start;
      front_circ *front = start;
      vertex_neighbor[i].insert(i);

      HV_circulator hv = IV[i]->vertex_begin();
      do {	
	front->next = new front_circ;
	front->next->prev = front; 
	front = front->next;
	front->hf = hv->prev()->opposite();

	//if (hv->facet()!=NULL) facet_occupied.insert(hv->facet()->index);
	vertex_neighbor[i].insert(hv->opposite()->vertex()->index);	

      } while (++hv != IV[i]->vertex_begin());
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
	if (facet_occupied.find(facet_index = front->hf->facet()->index) == facet_occupied.end() && !front->fr) {
	  // && vertex_outbound.find(vertex_index = front->hf->next()->vertex()->index) == vertex_outbound.end()) {

	  displace = IV[vertex_index]->point() - IV[i]->point();
	  if (displace * displace < square_distance_threshold) {
	    // insert new vertex to the front
	    facet_occupied.insert(facet_index);
	    vertex_neighbor[i].insert(vertex_index);
	    
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



      count += vertex_neighbor[i].size();
    }


    time(&end);
    std::cout << ", time: " << difftime( end, start) <<" seconds" << std::endl;
    
    return (double) count/(double) vertex_num;
  }






  int TriMesh::detect_vertex_keypoint(ScalarFunction &valueScalar, BooleanFunction &keyBoolean, int iter, int pre_iter){
    // in the implementation of keypoint detection, no threshold is taken.
    int count = 0;
    double coeff = .6;
    double neighbor_size = 3*coeff > 3.? 3*coeff: 3.;

    update_vertex_neighbor(neighbor_size); 

    // update static coefficient for smooth
    double sigma = coeff * avg_edge_len;
    std::vector<std::map<int, double> > coeff_list(vertex_num);
    for (int i=0; i<vertex_num; ++i) {
      Vector tmp;
      double scale, total_scale = 0;
      
      for (std::set<int>::iterator it = vertex_neighbor[i].begin();
	   it != vertex_neighbor[i].end(); ++it) {
	tmp = IV[*it]->point() - IV[i]->point();
	scale = CGAL::sqrt(tmp * tmp);
	scale = vertex_area[*it] * std::exp( - (scale * scale) / (2 * sigma * sigma));
	total_scale += scale;
	coeff_list[i][*it] = scale;
      }

      for (std::map<int, double>::iterator it = coeff_list[i].begin();
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
      //radio2 = 0.;

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
	    if (tmp)// && std::fabs(buffer_dv[one][j]) > radio2 )
	      { keyBoolean[j] = true; ++ count; continue;}

	  
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
	    if (tmp)// && std::fabs(buffer_dv[one][j])  > radio2 )
	      { keyBoolean[j] = true; ++ count; continue; }

	  }
	
	std::cout<< "#" << i << " detect Keypoint iteration\t Found " << count << std::endl;
	
      } else {
	std::cout <<"#" << i << " pre-smooth count" << std::endl;
      }      

      buffer_curr = (buffer_curr +1)%4;

    }

    valueScalar = buffer_v[ buffer_curr ];
    return count;
  }

}


