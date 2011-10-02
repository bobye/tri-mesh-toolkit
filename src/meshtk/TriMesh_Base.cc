/*
  FILE: TriMesh.cc This file is part of MeshTK.
  It is a C++ source file which implements the base functions of class TriMesh.
  
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

#include <iostream>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"

namespace meshtk {

  Curvature Curvature::tensor_compute(double e, double f, double g){
    double H= (e+g)/2;
    double G=std::sqrt(4*f*f+(e-g)*(e-g))/2;

    return Curvature(H+G, H-G, H, H*H-G*G);
  }

  TriMesh::TriMesh (){
    set_attribute_id = MESHTK_USER_ATTRIBUTE_START;
  };

  void TriMesh::read(std::string file, std::string type){

    if (type.compare("off")==0){
      std::ifstream mesh_Fin;
      file.append("."); file.append("off");    
      mesh_Fin.open(file.c_str());
      if (!mesh_Fin.is_open()){
	std::cerr << "Error: Cannot open file " << file << std::endl;
	exit(1);
      }

      mesh_Fin >> P;// mesh read
      mesh_Fin.close();
    }


    if (!P.is_pure_triangle()) {std::cerr<< "Error: input mesh has non-triangle facet!" <<std::endl; exit(1);}
    if (!P.is_closed()) {std::cout<< "Warning: input mesh seems not to be watertight" <<std::endl;} 
  };

  void TriMesh::write(std::string file, std::string type){

    if (type.compare("off")==0){
      std::ofstream mesh_Fout;
      file.append("."); file.append(type);
      mesh_Fout.open(file.c_str());
      mesh_Fout << P;// mesh output
      std::cout << "Export mesh to: " << file <<  std::endl;
      mesh_Fout.close();
    }
  };

  void TriMesh::init_index(){
    //int n=(P.size_of_halfedges()+P.size_of_border_edges())/2;

    halfedge_num = P.size_of_halfedges();
    vertex_num = P.size_of_vertices();
    facet_num = P.size_of_facets();

    IH = ISHalfedgeList(halfedge_num);
    IV =  ISVertexList(vertex_num);
    IF =  ISFacetList(facet_num);

    halfedge_vec = VectorFunction(halfedge_num);
    vertex_norm = VectorFunction(vertex_num);
    facet_norm = VectorFunction(facet_num);
    vertex_area = ScalarFunction(vertex_num);
    facet_area = ScalarFunction(facet_num);

    vertex_avg_len = ScalarFunction(vertex_num);

    int index_count=0;
    for(Vertex_iterator vitr= P.vertices_begin();vitr!= P.vertices_end();
	IV[index_count]=vitr, vitr->index = index_count++, vitr++);
    index_count=0;
    for(Halfedge_iterator eitr= P.halfedges_begin();eitr!= P.halfedges_end();
	IH[index_count]=eitr, eitr->index = index_count++, eitr++);
    index_count=0;
    for(Facet_iterator fitr= P.facets_begin(); fitr!= P.facets_end(); 
	IF[index_count]=fitr, fitr->index = index_count++, fitr++);
  };


  double TriMesh::update_halfedge(){

    double avg_len=0;
    Halfedge_handle h;
    for (int i=0;i<halfedge_num;i++) 
      {
	h = IH[i];
	halfedge_vec[i] = h->vertex()->point() - h->prev()->vertex()->point();
	avg_len += CGAL::sqrt(halfedge_vec[i] * halfedge_vec[i]);
      }
    return avg_edge_len = avg_len/halfedge_num;
  }


  double TriMesh::update_facet(){
    Vector normal;
    Halfedge_handle h;

    total_area=0;
    for (int i=0;i<facet_num;i++){
      h=IF[i]->halfedge();
      normal = CGAL::cross_product(halfedge_vec[h->index], halfedge_vec[h->next()->index]);
      total_area += (facet_area[i] = normal * normal) /2. ;
      facet_norm[i] = normal / CGAL::sqrt(facet_area[i]); 
      facet_area[i] /= 2.;
    }
    
    attribute[MESHTK_FACET_NORM] = &facet_norm;
  
    return  total_area;
  }


  void TriMesh::update_vertex(){
    /*
      double sigma = 2 * avg_edge_len;

      for (int i=0;i<vertex_num;i++){
      Vector normal(0,0,0), tmp;
      double scale;
      HV_circulator hv=IV[i]->vertex_begin();
      do{
      tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
      scale = CGAL::sqrt(tmp * tmp);
      normal = normal + std::exp(- scale * scale / (sigma *sigma)) * facet_norm[hv->facet()->index];
      }while (++hv!=IV[i]->vertex_begin());
      vertex_norm[i] = normal / CGAL::sqrt(normal * normal);
      }
    */


    for (int i=0;i<vertex_num;i++) {
      Point p= IV[i]->point();
      if (p.x() < coordinate_min_x || i==0) coordinate_min_x = p.x();
      if (p.y() < coordinate_min_y || i==0) coordinate_min_y = p.y();
      if (p.z() < coordinate_min_z || i==0) coordinate_min_z = p.z();
      if (p.x() > coordinate_max_x || i==0) coordinate_max_x = p.x();
      if (p.y() > coordinate_max_y || i==0) coordinate_max_y = p.y();
      if (p.z() > coordinate_max_z || i==0) coordinate_max_z = p.z();
    
    }

    for (int i=0; i<vertex_num; i++){
      HV_circulator hv=IV[i]->vertex_begin();
      double area = 0, total_len =0;
      Vector tmp;
      int k=0;

      do {
	if (hv->facet()==NULL) continue;
	area += facet_area[hv->facet()->index];

	tmp = halfedge_vec[hv->index];
	total_len += CGAL::sqrt(tmp * tmp);		    
	++k;      

      }while (++hv!=IV[i]->vertex_begin());
      vertex_area[i] = area/3.;
      vertex_avg_len[i] = total_len / k;

    }

    facet2vertex_point_average( facet_norm, vertex_norm, Vector(0,0,0));
    for (int i=0;i<vertex_num;i++){
      vertex_norm[i] = vertex_norm[i] / CGAL::sqrt(vertex_norm[i] * vertex_norm[i]);
    }
    
    attribute[MESHTK_VERTEX_NORM] = &vertex_norm;

    
  }


  void TriMesh::update_base(){//base update halfedge, facet, vertex.

    time_t  start, end; 
    time(&start);
    std::cout << "Update base ..."<< std::flush; 
    update_halfedge();

    update_facet();
  
    update_vertex();
    time(&end);
    std::cout << "\t time: " << difftime( end, start) <<" seconds" << std::endl;

  };
  

  TriMesh::~TriMesh(){    
  }


}