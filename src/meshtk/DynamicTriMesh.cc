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
#include <algorithm>


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
    for (int i = 0; i < vertex_num; i++){
      Vector tmp, vec(0,0,0);
      double scale, total_scale = vertex_area[i], sigma = coeff * avg_edge_len /2.;
      HV_circulator hv=IV[i]->vertex_begin();
            

      do {
	//	sigma = coeff * vertex_avg_len[i];

	if (hv->facet() == NULL) continue;
	
	tmp = ( - halfedge_vec[hv->index] + halfedge_vec[hv->next()->index] )/3.;
	//tmp2 = tmp - (tmp * vertex_norm[i]) * vertex_norm[i];
	//scale = CGAL::sqrt(tmp2 * tmp2);

	scale = CGAL::sqrt(tmp * tmp);

	scale = facet_area[hv->facet()->index] * std::exp (- (scale * scale) / (2 * sigma * sigma));
	total_scale +=scale;

	vec = vec + scale * tmp;
 
      }while (++hv != IV[i]->vertex_begin());
      
      vertex_coord[i] = IV[i]->point() + vec / total_scale;
      
    }

    restore_coord();
  }


  void DynamicTriMesh::update_vertex_salient(int iter, int pre_iter){
    vertex_salient.resize(vertex_num);
    vertex_salient_sup.resize(vertex_num);
    vertex_salient_inf.resize(vertex_num);

    // buffers of mean curvature
    ScalarFunction buffer_hcurv[4]; 
    ScalarFunction buffer_doh[4];

    for (int i=0;i<4;i++ ) {buffer_hcurv[i].resize(vertex_num); buffer_doh[i].resize(vertex_num);}

    double coeff = 1.;
    int buffer_curr=0, buffer_renew=0;
    
    //pre-smooth, if the input mesh is manifold mesh, set pre_iter as default
    while (pre_iter-- > 0)  { gaussian_smooth(1.); update_base();}


    update_curvature();    
    buffer_hcurv[buffer_renew++] = vertex_hcurv;

    gaussian_smooth(coeff);
    update_base();
    update_curvature();
    buffer_hcurv[buffer_renew++] = vertex_hcurv;

    gaussian_smooth(coeff);
    update_base();
    update_curvature();
    buffer_hcurv[buffer_renew++] = vertex_hcurv;

    for (int i=0; i < vertex_num; ++i) {
      buffer_doh[0][i] = buffer_hcurv[1][i] - buffer_hcurv[0][i];
      buffer_doh[1][i] = buffer_hcurv[2][i] - buffer_hcurv[1][i];
    }

    for (int i=0; i < iter; ++i) {
      
      gaussian_smooth(coeff);
      update_base();
      update_curvature();
      buffer_hcurv[buffer_renew] = vertex_hcurv;
      
      int zero = buffer_curr, 
	one = (buffer_curr +1 )%4,
	two = (buffer_curr +2 )%4,
	three = (buffer_curr +3)%4;
      
      for (int j=0; j < vertex_num; ++j) buffer_doh[two][j] = buffer_hcurv[three][j] - buffer_hcurv[two][j];
            
      for (int j=0; j < vertex_num; ++j) {// the algorithm below is not much efficient
	
	HV_circulator hv=IV[j]->vertex_begin();	

	bool tmp = true;
	if (buffer_doh[one][j] < buffer_doh[two][j] || buffer_doh[one][j] < buffer_doh[zero][j])
	  tmp = false;
	
	do {
	  if (buffer_doh[one][j] < buffer_doh[one][hv->next()->vertex()->index] ||
	      buffer_doh[two][j] < buffer_doh[two][hv->next()->vertex()->index] ||
	      buffer_doh[zero][j] < buffer_doh[zero][hv->next()->vertex()->index]) {
	    tmp = false;
	    break;
	  }	  
	} while (++hv != IV[j]->vertex_begin());
	vertex_salient_sup[j] = vertex_salient_sup[j] || tmp;
	
	tmp = true;
	if (buffer_doh[one][j] > buffer_doh[two][j] || buffer_doh[one][j] > buffer_doh[zero][j])
	  tmp = false;
	
	do {
	  if (buffer_doh[one][j] > buffer_doh[one][hv->next()->vertex()->index] ||
	      buffer_doh[two][j] > buffer_doh[two][hv->next()->vertex()->index] ||
	      buffer_doh[zero][j] > buffer_doh[zero][hv->next()->vertex()->index]) {
	    tmp = false;
	    break;
	  }	  
	} while (++hv != IV[j]->vertex_begin());
	vertex_salient_inf[j] = vertex_salient_inf[j] || tmp;

	vertex_salient[j] = vertex_salient[j] || vertex_salient_sup[j] || vertex_salient_inf[j];

      }
      std::cout<< "Salient Detection Iteration: "<< i << "\t Found " << std::accumulate( vertex_salient.begin(), vertex_salient.end(), 0) <<" salient points" << std::endl;

      buffer_renew = (buffer_renew +1) % 4;
      buffer_curr = (buffer_curr +1) % 4;
    }
    
  }

}
