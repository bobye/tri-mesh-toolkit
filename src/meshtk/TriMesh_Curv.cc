/*
  FILE: TriMesh_Curv.cc This file is part of MeshTK.
  It is a C++ source file which implements the curvature computation 
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

namespace meshtk {

  void TriMesh::update_vertex_localchart(){
    vertex_LC[0] = VectorFunction(vertex_num);
    vertex_LC[1] = VectorFunction(vertex_num);
  
    for (int i=0;i<vertex_num;i++) localchart(vertex_LC[0][i],vertex_LC[1][i], vertex_norm[i]);
  }

  void TriMesh::update_facet_localchart(){
    facet_LC[0] = VectorFunction(facet_num);
    facet_LC[1] = VectorFunction(facet_num);

    for (int i=0;i<facet_num;i++) localchart(facet_LC[0][i],facet_LC[1][i], facet_norm[i]);
  }



  void TriMesh::update_facet_curvature(){

    facet_CT[0].resize(facet_num);
    facet_CT[1].resize(facet_num);
    facet_CT[2].resize(facet_num);

    facet_PC0.resize(facet_num);
    facet_PC1.resize(facet_num);
    facet_hcurv.resize(facet_num);
    facet_kcurv.resize(facet_num);
  
    //for (int i=0;i<1;i++){
    for (int i=0;i<facet_num;i++){

      //compute curvature for each face

      //step 1. assembly matrix
      Halfedge_handle h[3];
      h[0]=IF[i]->halfedge();
      h[1]=h[0]->next();
      h[2]=h[1]->next();

      Vector e[3],n[3];
      for (int j=0;j<3;j++) e[j] = halfedge_vec[h[j]->index];
    
      for (int j=0;j<3;j++) n[j] = vertex_norm[h[j]->vertex()->index] - vertex_norm[h[j]->prev()->vertex()->index];

      double e_LC[3][2], n_LC[3][2];
      for (int j=0;j<3;j++) {
	localcoord(e[j],facet_LC[0][i],facet_LC[1][i], e_LC[j]);
	localcoord(n[j],facet_LC[0][i],facet_LC[1][i], n_LC[j]);
      }


      double matrix_l[2][2]={}, matrix_r[2][2]={};
      for (int j=0;j<2;j++)
	for (int k=0;k<2;k++)
	  for (int l=0;l<3;l++){
	    matrix_l[j][k]+=e_LC[l][j]*e_LC[l][k];
	    matrix_r[j][k]+=n_LC[l][j]*e_LC[l][k];
	  }
    
      double final_linv[3][3]={}, final_r[3];
      double a,b,c, deter;
      a = matrix_l[0][0]; b = matrix_l[0][1]; c = matrix_l[1][1];
      deter=(a+c)*(a*c-b*b);
    
      final_r[0] = matrix_r[0][0];
      final_r[1] = matrix_r[0][1]+matrix_r[1][0];
      final_r[2] = matrix_r[1][1];

      //step 2. solve matrix    
      final_linv[0][0]=(a+c)*c-b*b;
      final_linv[0][1]=final_linv[1][0]=-b*c;
      final_linv[0][2]=final_linv[2][0]=b*b;
      final_linv[1][1]=a*c;
      final_linv[1][2]=final_linv[2][1]=-b*a;
      final_linv[2][2]=(a+c)*a-b*b;

    
      facet_CT[0][i] = (final_linv[0][0]*final_r[0] + final_linv[0][1]*final_r[1] + final_linv[0][2]*final_r[2])/deter;
    
      facet_CT[1][i] = (final_linv[1][0]*final_r[0] + final_linv[1][1]*final_r[1] + final_linv[1][2]*final_r[2])/deter;
      facet_CT[2][i] = (final_linv[2][0]*final_r[0] + final_linv[2][1]*final_r[1] + final_linv[2][2]*final_r[2])/deter;        


      Curvature curv = Curvature::tensor_compute(facet_CT[0][i], facet_CT[1][i], facet_CT[2][i]);

      facet_PC0[i] = curv.principle_curv0;
      facet_PC1[i] = curv.principle_curv1;
      facet_hcurv[i] = curv.mean_curv;
      facet_kcurv[i] = curv.gaussian_curv;

    }

    attribute[MESHTK_FACET_PC0] = &facet_PC0;
    attribute[MESHTK_FACET_PC1] = &facet_PC1;
    attribute[MESHTK_FACET_HCURV] = &facet_hcurv;
    attribute[MESHTK_FACET_KCURV] = &facet_kcurv;    

  }



  void TriMesh::update_vertex_curvature(){
    vertex_CT[0].resize(vertex_num);
    vertex_CT[1].resize(vertex_num);
    vertex_CT[2].resize(vertex_num);

    vertex_PC0.resize(vertex_num);
    vertex_PC1.resize(vertex_num);

    vertex_hcurv.resize(vertex_num);
    vertex_kcurv.resize(vertex_num);
    double sigma = avg_edge_len;
    Vector tmp;

    for (int i=0; i < vertex_num; ++i) {
      HV_circulator hv = IV[i]->vertex_begin();
      double total_scale=0, scale;
      Vector tmp;
      vertex_CT[0][i] = vertex_CT[1][i] = vertex_CT[2][i] =0;

      sigma = vertex_avg_len[i];

      do{
	if (hv->facet() == NULL) continue;

	double x[2][2]; 
	int j=hv->facet()->index;

	tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
	scale = CGAL::sqrt(tmp * tmp) /3.;

	scale = facet_area[j]*std::exp(- scale * scale / (2 * sigma * sigma));
	total_scale += scale;


	for (int m=0;m<2;++m) for (int k=0;k<2;++k) x[m][k] = vertex_LC[m][i] * facet_LC[k][j];

      
	vertex_CT[0][i] += scale * 
	  (facet_CT[0][j]*x[0][0]*x[0][0] + 2*facet_CT[1][j]*x[0][0]*x[0][1] + facet_CT[2][j]*x[0][1]*x[0][1]);
	vertex_CT[1][i] += scale *
	  (facet_CT[0][j]*x[0][0]*x[1][0] + facet_CT[1][j]*(x[0][0]*x[1][1] + x[0][1]*x[1][0])+ facet_CT[2][j]*x[0][1]*x[1][1]);
	vertex_CT[2][i] += scale *
	  (facet_CT[0][j]*x[1][0]*x[1][0] + 2*facet_CT[1][j]*x[1][0]*x[1][1] + facet_CT[2][j]*x[1][1]*x[1][1]);

      }while (++hv != IV[i]->vertex_begin());

      vertex_CT[0][i]/=total_scale;
      vertex_CT[1][i]/=total_scale;
      vertex_CT[2][i]/=total_scale;

 
      Curvature curv = Curvature::tensor_compute(vertex_CT[0][i], vertex_CT[1][i], vertex_CT[2][i]);

      vertex_PC0[i] = curv.principle_curv0;
      vertex_PC1[i] = curv.principle_curv1;
      vertex_hcurv[i] = curv.mean_curv;
      vertex_kcurv[i] = curv.gaussian_curv;
    }

    attribute[MESHTK_VERTEX_PC0] = &vertex_PC0;
    attribute[MESHTK_VERTEX_PC1] = &vertex_PC1;
    attribute[MESHTK_VERTEX_HCURV] = &vertex_hcurv;
    attribute[MESHTK_VERTEX_KCURV] = &vertex_kcurv;
    /*
      facet2vertex_point_average( facet_hcurv, vertex_hcurv, 0.);
      facet2vertex_point_average( facet_kcurv, vertex_kcurv, 0.);
    */
  }

  void TriMesh::update_curvature(){

    time_t  start, end; 
    time(&start);
    std::cout << "Update curvature ..."<< std::flush; 

    update_facet_localchart();
    update_facet_curvature();

    update_vertex_localchart();
    update_vertex_curvature();

    time(&end);
    std::cout << "\t[done] " << difftime( end, start) <<" seconds" << std::endl;


  }

}
