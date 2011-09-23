/*
  FILE: meshtk/TriMesh.hh This file is part of MeshTK.
  It is a C++ header file which defines the interface of class TriMesh. 
  
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

#ifndef _TRIMESH_HH_
#define _TRIMESH_HH_
#include <vector>
#include "mesh_topo.hh"

namespace meshtk {

  typedef std::vector<Halfedge_handle> ISHalfedge_list;
  typedef std::vector<Vertex_handle> ISVertex_list;
  typedef std::vector<Facet_handle> ISFacet_list;

  typedef std::vector<Vector> Vec_Fun;
  typedef std::vector<double> Scalar_Fun;
  typedef std::vector<double> Bool_Fun;// For marker.


  struct Curvature {
    Curvature (const double pc0, const double pc1, const double hc, const double kc) :PC0(pc0), PC1(pc1), hcurv(hc), kcurv(kc){}
    const double PC0, PC1, hcurv, kcurv;
  };

  class TriMesh {//topo ref system //index system of items of polyhedron

    Polyhedron P;

    double coord_min_x, coord_min_y, coord_min_z;
    double coord_max_x, coord_max_y, coord_max_z;


    Vec_Fun vertex_norm;//normal vector
    Vec_Fun facet_norm;


    ISHalfedge_list IH; // 
    ISVertex_list IV;
    ISFacet_list IF;

    //item attributes
    Vec_Fun halfedge_vec;    


    Vec_Fun vertex_LC[2];//LC: local chart
    Vec_Fun facet_LC[2];
  
    Scalar_Fun vertex_CT[3];//CT: cuvature tensor [ 0 1 ]
    Scalar_Fun facet_CT[3];//                     [ 1 2 ]

    Scalar_Fun vertex_avg_len;
    
    Scalar_Fun vertex_area;
    Scalar_Fun facet_area;

    Scalar_Fun vertex_PC0;//PC: principle cuvature
    Scalar_Fun vertex_PC1;
    Scalar_Fun facet_PC0;
    Scalar_Fun facet_PC1;

    Scalar_Fun vertex_hcurv;//mean curvature
    Scalar_Fun vertex_kcurv;//Gaussian curvature

    Scalar_Fun facet_hcurv;//mean curvature
    Scalar_Fun facet_kcurv;//Gaussian curvature


    double avg_edge_len;


    double update_halfedge();
    double update_facet();
    void update_vertex();

    void update_facet_localchart();
    void update_vertex_localchart();


    void update_facet_curvature();
    void update_vertex_curvature();

    std::map<unsigned, Scalar_Fun *> feature_scalar_fun;
    std::map<unsigned, Vec_Fun *> feature_vec_fun;

  public:  

    double total_area;
    int halfedge_num;
    int vertex_num;
    int facet_num;

    static const unsigned VERTEX_PC0 = 0;
    static const unsigned VERTEX_PC1 = 1;
    static const unsigned VERTEX_HCURV = 2;
    static const unsigned VERTEX_KCURV = 3;

    static const unsigned FACET_PC0 = 4;
    static const unsigned FACET_PC1 = 5;
    static const unsigned FACET_HCURV = 6;
    static const unsigned FACET_KCURV = 7;


    void read(std::string, std::string);
    void write(std::string, std::string);
    void init_index();
    void update_base();
    void update_curvature();

    Scalar_Fun * feature_scalar_extract(unsigned );
    Vec_Fun * feature_vec_extract(unsigned );

    template <class T>
    void facet2vertex_area_average(std::vector<T> &f, std::vector<T> &v, T zero){
    };

    template <class T>
    void facet2vertex_point_average(std::vector<T> &f, std::vector<T> &v, T zero){
      double sigma = avg_edge_len;

      for (int i=0;i<vertex_num;i++){
	Vector tmp;
	double scale, total_scale=0;
	HV_circulator hv=IV[i]->vertex_begin();
	T vec = zero;

	sigma = vertex_avg_len[i] *2. /3.;
	
	do{
	  if (hv->facet()==NULL) continue;

	  tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
	  scale = CGAL::sqrt(tmp * tmp) /3.;		    
	  scale = facet_area[hv->facet()->index] * std::exp(- scale * scale / (2 * sigma * sigma));	    

	  total_scale += scale;
	  vec = vec + scale * f[hv->facet()->index];

	} while (++ hv!=IV[i]->vertex_begin());


	v[i] = vec/total_scale;
      }
    };

    template <class T>
    void vertex2facet_average(std::vector<T> &v, std::vector<T> &f){
      for (int i=0;i<facet_num;i++){
	// ...
      }
    };






    friend class MeshPainter;
  

  };

}
#endif /* _TRIMESH_HH_ */
