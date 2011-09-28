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

namespace meshtk{


  double TriMesh::update_vertex_neighbor(double coeff){
    double distance_threshold = coeff * avg_edge_len; 
    vertex_neighbor.resize(vertex_num);
    int count=0;

    // The algorithm implemented is a wide search 
    
    for (int i = 0; i < vertex_num; ++i) {
      
      std::set<int> buffer;
      buffer.insert(i);

      while (!buffer.empty()) {
	std::set<int>::iterator it = buffer.begin();
	Vector displace = IV[*it]->point()- IV[i]->point();
	
	if (displace * displace < distance_threshold * distance_threshold) {
	  vertex_neighbor[i].insert(*it);
	  vertex_neighbor[*it].insert(i);
	  
	  HV_circulator hv = IV[*it]->vertex_begin();
	  do {
	    int j =hv->opposite()->vertex()->index;
	    if (vertex_neighbor[i].count(j) == 0) buffer.insert(j);
	  }while (++hv != IV[*it]->vertex_begin());
	}

	buffer.erase(it);
	
      }

      count += vertex_neighbor[i].size();

    }

    return (double )count / (double) vertex_num;
  }


  void TriMesh::detect_vertex_keypoint(ScalarFunction &valueScalar, BooleanFunction &keyBoolean, int iter, int pre_iter){
    double coeff = 1.4;
    
    update_vertex_neighbor(3 * coeff);

    ScalarFunction buffer_v[4];
    ScalarFunction buffer_dv[4];

    for (int i=0; i<4; ++i) {buffer_v[i].resize(vertex_num); buffer_dv[i].resize(vertex_num);}
    int buffer_curr = 0;
    
    buffer_v[0] = valueScalar;

    gaussian_smooth_vertex(coeff, buffer_v[0], buffer_v[1], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[1], buffer_v[2], 0.);
    gaussian_smooth_vertex(coeff, buffer_v[2], buffer_v[3], 0.);

    for (int i = 0; i < vertex_num; ++i) {
      buffer_dv[0][i] = buffer_v[1][i] - buffer_v[0][i];
      buffer_dv[1][i] = buffer_v[2][i] - buffer_v[1][i];
    }
    
    for (int i=0; i< iter; ++i) {
      
      int zero = buffer_curr, one=(buffer_curr +1)%4, 
	two=(buffer_curr +2)%4, three=(buffer_curr +3)%4;

      gaussian_smooth_vertex(coeff, buffer_v[two], buffer_v[three], 0.);
      for (int j=0; j < vertex_num; ++j) buffer_dv[two][j] = buffer_v[three][j] - buffer_v[two][j];

      if (pre_iter--<=0) {
	for (int j=0; j < vertex_num; ++j) {
	  HV_circulator hv=IV[j]->vertex_begin();
	  bool tmp = true;
	  
	  if (buffer_dv[one][j] < buffer_dv[two][j] || buffer_dv[one][j] < buffer_dv[zero][j])
	    tmp = false;
	  
	  do {
	    if (buffer_dv[one][j] < buffer_dv[one][hv->next()->vertex()->index] ||
		buffer_dv[two][j] < buffer_dv[two][hv->next()->vertex()->index] ||
		buffer_dv[zero][j] < buffer_dv[zero][hv->next()->vertex()->index]) {
	      tmp = false; break;
	    }
	  }while (++hv!=IV[j]->vertex_begin());
	  
	  keyBoolean[j] = keyBoolean[j] || tmp;
	  
	  hv=IV[j]->vertex_begin(); tmp = true;
	  if (buffer_dv[one][j] > buffer_dv[two][j] || buffer_dv[one][j] > buffer_dv[zero][j])
	    tmp = false;
	  
	  do {
	    if (buffer_dv[one][j] > buffer_dv[one][hv->next()->vertex()->index] ||
		buffer_dv[one][j] > buffer_dv[two][hv->next()->vertex()->index] ||
		buffer_dv[one][j] > buffer_dv[zero][hv->next()->vertex()->index]) {
	      tmp = false; break;
	    }
	  }while (++hv!=IV[j]->vertex_begin());
	  
	  keyBoolean[j] = keyBoolean[j] || tmp;	  
	}
	
	std::cout<< "Keypoint detection iteration: "<<i<< "\t Found " << std::accumulate(keyBoolean.begin(), keyBoolean.end(), 0)<< std::endl;
	
      } else {
	std::cout <<"Smooth count"<<std::endl;
      }      

      buffer_curr = (buffer_curr +1)%4;

    }
    
    
  }

}
