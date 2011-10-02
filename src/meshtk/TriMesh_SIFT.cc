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

    double square_distance_threshold = coeff * avg_edge_len; 
    square_distance_threshold *=square_distance_threshold;

    vertex_neighbor.clear();
    vertex_neighbor.resize(vertex_num);
    int count=0;// to record average neighbor size of vertices

    // The algorithm implemented is a wide search 
    // which is not very efficient
    /*
    for (int i = 0; i < vertex_num; ++i) {
      int j;
      std::set<int>::iterator it;
      std::set<int> buffer;
      Vector displace;
      buffer.insert(i);
      HV_circulator hv;

      while (!buffer.empty()) {
	it = buffer.begin();
	displace = IV[*it]->point()- IV[i]->point();
	
	if (displace * displace < square_distance_threshold) {
	  vertex_neighbor[i].insert(*it);
	  
	  hv = IV[*it]->vertex_begin();
	  do {
	    j =hv->opposite()->vertex()->index;
	    if (vertex_neighbor[i].find(j) == vertex_neighbor[i].end()) buffer.insert(j);
	  }while (++hv != IV[*it]->vertex_begin());
	}

	buffer.erase(it);	
      }

      count += vertex_neighbor[i].size();
    }
    */

    // the algorithm implemented is a front propagation one
    // initialization from a 1-ring neighbor
    struct front_circ {
      Halfedge_handle hf;
      front_circ *prev;
      front_circ *next;
    };
    

    for (int i = 0; i < vertex_num; ++i){      
      Vector displace;
      front_circ *start = new front_circ; start->prev = start->next = start;
      front_circ *front;
      vertex_neighbor[i].insert(i);

      HV_circulator hv = IV[i]->vertex_begin();
      int M=-1;

      do {	
	front = new front_circ;
	front->hf = hv->prev()->opposite();
	front->next = start; front->prev = start->prev;
	start->prev->next = front; start->prev = front; 
	vertex_neighbor[i].insert(hv->opposite()->vertex()->index);
	if (i==M) std::cout<< hv->opposite()->vertex()->index <<" -> "<< hv->prev()->opposite()->vertex()->index << std::endl;
	//getchar();
	
      } while (++hv != IV[i]->vertex_begin());
      start->prev->next = start->next; start->next->prev = start->prev;
      delete start;
      start = front;
      bool front_reached = false;
      
      do {	  
	if (i==M) getchar();
	
	front_reached = true;
	//while (front->hf->opposite() == front->prev->hf || front->hf->opposite() == front->next->hf ){
	while (front->hf->opposite() == front->next->hf) {	  
	  // delete verbose pair of halfedges to improve performance and 
	  // insert the corresponding vertex into neighbor set
	  if (i==M) std::cout<<"- "<<front->hf->vertex()->index <<std::endl;

	  front->next = front->next->next; 
	  delete front->next->prev;
	  front->next->prev = front; 
	  front = front->prev;
	  front->next = front->next->next; 
	  delete front->next->prev;
	  front->next->prev = front;
	  start = front;
	}

	/*
	while (front->hf->opposite() == front->prev->hf) {
	  front = front->prev;
	  front->next = front->next->next; 
	  delete front->next->prev;
	  front->next->prev = front; 
	  front = front->prev;
	  front->next = front->next->next; 
	  delete front->next->prev;
	  front->next->prev = front;	  
	  start = front;	  
	}
	//	}
	*/
	if (front->hf->facet()==NULL){
	  // delete verbose boundary halfedge to improve performance and 
	  // insert the corresponding vertex into neighbor set
	  //std::cout << "d " << front->hf->opposite()->vertex()->index << std::endl;
	  front = front->prev;
	  front->next = front->next->next;
	  delete front->next->prev;
	  front->next->prev = front;
	  start = front;
	}

	// to test whether propagate toward the facet
	displace = front->hf->next()->vertex()->point() - IV[i]->point();
	if (displace * displace < square_distance_threshold) {
	  // insert new vertex to the front
	  
	  vertex_neighbor[i].insert(front->hf->next()->vertex()->index);
	  front_circ *insert = new front_circ;
	  insert->hf = front->hf->prev()->opposite();
	  insert->prev = front->prev; insert->next = front;
	  front->hf = front->hf->next()->opposite();
	  front->prev->next = insert; front->prev = insert;
	  front_reached =false;
	  start=front;
	  if (i==M) std::cout<<"+ " <<front->hf->opposite()->vertex()->index << " -> " <<front->hf->vertex()->index << std::endl;
	}
	else {
	  front = front->next; 
	  if (i==M) std::cout<<"h "<<front->hf->vertex()->index <<std::endl;
	}
	  

      } while (front!= start || !front_reached);
      
      front->next->prev = NULL;
      while (front->prev!= NULL) { front = front->prev; delete front->next;  } delete front;

      std::cout << i << std::endl;

      count += vertex_neighbor[i].size();
    }
    
    return (double) count/(double) vertex_num;
  }






  int TriMesh::detect_vertex_keypoint(ScalarFunction &valueScalar, BooleanFunction &keyBoolean, int iter, int pre_iter){
    // in the implementation of keypoint detection, no threshold is taken.
    int count = 0;
    double coeff = .6;
    double neighbor_size = 3*coeff > 3.? 3*coeff: 3.;
    
    time_t  start, end; time(&start);
    std::cout << "Update neighbors, size: " << neighbor_size << std::flush; 
    update_vertex_neighbor(neighbor_size); time(&end);
    std::cout << ", time: " << difftime( end, start) <<" seconds" << std::endl;

    ScalarFunction buffer_v[4];
    ScalarFunction buffer_dv[4];

    for (int i=0; i<4; ++i) {buffer_v[i].resize(vertex_num); buffer_dv[i].resize(vertex_num);}
    int buffer_curr = 0;
    
    buffer_v[0] = valueScalar;

    gaussian_smooth_vertex(coeff, buffer_v[0], buffer_v[1], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[1], buffer_v[2], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[2], buffer_v[3], 0.);

    for (int i = 0; i < vertex_num; ++i) {
      buffer_dv[0][i] =   (buffer_v[1][i] - buffer_v[0][i]) / std::log(2.);
      buffer_dv[1][i] =   (buffer_v[2][i] - buffer_v[1][i]) / std::log(1.5);
    }
    
    for (int i=0; i< iter; ++i) {
      
      //      double radio = 0,// (0.001 / std::sqrt(total_area)) * coeff * coeff, 
      //radio2 = 0.;

      int zero = buffer_curr, one=(buffer_curr +1)%4, 
	two=(buffer_curr +2)%4, three=(buffer_curr +3)%4;

      gaussian_smooth_vertex(coeff, buffer_v[two], buffer_v[three], 0.);
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
	
	std::cout<< "Keypoint detection iteration: "<<i<< "\t Found " << count << std::endl;
	
      } else {
	std::cout <<"Smooth count"<<std::endl;
      }      

      buffer_curr = (buffer_curr +1)%4;

    }

    valueScalar = buffer_v[ buffer_curr ];
    return count;
  }

}


