/*
  FILE: DynamicTriMesh.hh This file is part of tri-mesh-toolkit.
  It is a C++ source file which implemented class DynamicTriMesh, achieving smoothing, 
  denoising and remeshing.
  
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

#include "meshtk/DynamicTriMesh.hh"
#include "meshtk/mesh_assist.hh"
#include <algorithm>


namespace meshtk {
  
  void DynamicTriMesh::init_index(){
    TriMesh::init_index();
    
    vertex_coord = PointFunction(vertex_num);
    facet_coord = PointFunction(facet_num);
  }

  void DynamicTriMesh::load_coord(){    
    vertex_coord.resize(vertex_num);


    for (int i=0;i<vertex_num; i++){
      vertex_coord[i] =  IV[i]->point();
    }

    attribute[MESHTK_VERTEX_COORD] = &vertex_coord;

  }
  
  void DynamicTriMesh::restore_coord(){
    for (int i=0;i<vertex_num; i++){
      IV[i]->point() = vertex_coord[i];
    }
  }


  void DynamicTriMesh::add_mesh_noise(double coeff) {

    double move = coeff * avg_edge_len;

    load_coord();
    srand((unsigned)time(NULL));
    for (int i = 0;i<vertex_num; i++) 
      vertex_coord[i] = vertex_coord[i] + ((double )rand()/(double) RAND_MAX -.5) * move * vertex_norm[i];
    restore_coord();

  }

  void DynamicTriMesh::gaussian_smooth(double coeff){
    // preconditioned with vertex_neighbor

    //update_vertex_neighbor(3 * coeff);

    VectorFunction smoothed_normal(vertex_num);

    update_vertex_neighbor_euclidean(3 * coeff);
    neighbor_distance_map  = & vertex_neighbor_euclidean;

    gaussian_smooth_vertex(coeff, vertex_norm, smoothed_normal, Vector(0,0,0));
    for (int j=0; j<vertex_num; j++) vertex_norm[j] = smoothed_normal[j]/CGAL::sqrt(smoothed_normal[j] * smoothed_normal[j]);
      
    update_curvature();



    for (int i = 0; i < vertex_num; i++){
      Vector tmp;
      double tmp_scalar, norm_displace=0;

      double scale, total_scale = vertex_area[i];

      double sigma = vertex_avg_len[i];
      double coord[2];

      HV_circulator hv = IV[i]->vertex_begin();
      do {

	tmp = hv->opposite()->vertex()->point() - IV[i]->point(); 

	localcoord(tmp, vertex_LC[0][i], vertex_LC[1][i], coord);
	tmp_scalar = tmp * vertex_norm[i] + (vertex_CT[0][i]*coord[0]*coord[0] + 2*vertex_CT[1][i]*coord[0]*coord[1] + vertex_CT[2][i]*coord[1]*coord[1]);
	
	scale = CGAL::sqrt(tmp * tmp);
	scale = vertex_area[hv->opposite()->vertex()->index] * std::exp( - (scale * scale) / (2 * sigma * sigma));
	total_scale += scale;
	norm_displace += scale * tmp_scalar;


      }while (++hv!=IV[i]->vertex_begin());

      vertex_coord[i] = IV[i]->point() + (norm_displace / total_scale) * vertex_norm[i];
      
    }
      

    restore_coord();
  }


  void DynamicTriMesh::detect_vertex_salient(int iter, int pre_iter){
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
    while (pre_iter-- > 0)  { 
      gaussian_smooth(coeff); update_base();
    }


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
      buffer_doh[0][i] =  (buffer_hcurv[1][i] - buffer_hcurv[0][i]);
      buffer_doh[1][i] =  (buffer_hcurv[2][i] - buffer_hcurv[1][i]);
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
      
      for (int j=0; j < vertex_num; ++j) 
	buffer_doh[two][j] =  (buffer_hcurv[three][j] - buffer_hcurv[two][j]);
            
      for (int j=0; j < vertex_num; ++j) {// the algorithm below is not much efficient
	
	/****************************************************************************/
	HV_circulator hv=IV[j]->vertex_begin();	
	bool tmp = true;
	if (buffer_doh[one][j] < buffer_doh[two][j] || buffer_doh[one][j] < buffer_doh[zero][j])
	  tmp = false;	
	else do {
	    if (buffer_doh[one][j] < buffer_doh[one][hv->next()->vertex()->index] ||
		buffer_doh[two][j] < buffer_doh[two][hv->next()->vertex()->index] ||
		buffer_doh[zero][j] < buffer_doh[zero][hv->next()->vertex()->index]) {
	      tmp = false;
	      break;
	    }	  
	  } while (++hv != IV[j]->vertex_begin());
	
	
	vertex_salient_sup[j] = vertex_salient_sup[j] || tmp;

	/****************************************************************************/
	hv = IV[j]->vertex_begin();
	tmp = true;
	if (buffer_doh[one][j] > buffer_doh[two][j] || buffer_doh[one][j] > buffer_doh[zero][j])
	  tmp = false;	
	else do {
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


      std::cout<< "Salient Detection Iteration: "<< i << "\t Found " << std::accumulate( vertex_salient.begin(), vertex_salient.end(), 0) <<
	"("<< std::accumulate( vertex_salient_sup.begin(), vertex_salient_sup.end(), 0) <<
	","<< std::accumulate( vertex_salient_inf.begin(), vertex_salient_inf.end(), 0) << 
	") salient points" << std::endl;

      buffer_renew = (buffer_renew +1) % 4;
      buffer_curr = (buffer_curr +1) % 4;
    }
   
    attribute[MESHTK_VERTEX_SALIENT] = &vertex_salient;
    attribute[MESHTK_VERTEX_SALIENT_SUP] = &vertex_salient_sup;
    attribute[MESHTK_VERTEX_SALIENT_INF] = &vertex_salient_inf;
    
  }


  void DynamicTriMesh::remove_mesh_facets(int count, double percentage) {
    int pcount = (int) facet_num * percentage;
    if (count<0 || count > pcount) count = pcount;
    
    for (int i = 0; i < count; ++i) {
      Halfedge_iterator h = P.halfedges_begin();
      int itercount = rand() % (halfedge_num/2);
      while (itercount>0) {++h;--itercount;}
      while (h->is_border()) ++h;
      P.erase_facet(h);
    }
      
    P.keep_largest_connected_components(1);
  }



}


