//#include "mesh_topo.h"
#include "TriMesh.h"
#include <iostream>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>


#include "mesh_assist.h"

void TriMesh::read(std::string file, std::string type){
  if (type.compare("off")==0){
    std::ifstream mesh_Fin;
    file.append("."); file.append("off");    
    mesh_Fin.open(file.c_str());
    if (!mesh_Fin.is_open())
      std::cerr << "Cannot open file " << file << std::endl;
    mesh_Fin >> P;// mesh read
    mesh_Fin.close();
  }
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
  int n=(P.size_of_halfedges()+P.size_of_border_edges())/2;

  edge_num = n;
  vertex_num=P.size_of_vertices();
  facet_num=P.size_of_facets();

  IH = ISHalfedge_list(2*n);
  IV =  ISVertex_list(vertex_num);
  IF =  ISFacet_list(facet_num);

  halfedge_vec.resize(2*n);
  facet_norm.resize(P.size_of_facets());
  vertex_norm.resize(P.size_of_vertices());
  facet_area.resize(P.size_of_facets());

  int index_count=0;
  for(Vertex_iterator vitr= P.vertices_begin();vitr!= P.vertices_end();
      IV[index_count]=vitr, vitr->index = index_count++, vitr++);
  index_count=0;
  for(Edge_iterator eitr= P.edges_begin();eitr!= P.edges_end();
      IH[index_count]=eitr, IH[index_count+n]=eitr->opposite(),
	eitr->index = index_count, eitr->opposite()->index = index_count +n,
	index_count++, eitr++);
  index_count=0;
  for(Facet_iterator fitr= P.facets_begin(); fitr!= P.facets_end(); IF[index_count]=fitr, fitr->index=index_count++,fitr++);
};


double TriMesh::update_halfedge(){

  double avg_len=0;
  Halfedge_handle h;
  for (int i=0;i<2*edge_num;i++) 
    if ((h=IH[i])!=NULL){    
      halfedge_vec[i] = h->vertex()->point() - h->opposite()->vertex()->point();
      avg_len += CGAL::sqrt(halfedge_vec[i] * halfedge_vec[i]);
    }
  return avg_edge_len = avg_len/(2*edge_num-P.size_of_border_edges());
}


double TriMesh::update_facet(){
  Vector normal;
  Halfedge_handle h;

  total_area=0;
  for (int i=0;i<facet_num;i++){
    h=IF[i]->halfedge();
    normal = CGAL::cross_product(halfedge_vec[h->index], halfedge_vec[h->next()->index]);
    total_area += (facet_area[i] = normal * normal)/2. ;
    facet_norm[i] = normal / CGAL::sqrt(facet_area[i]); 
    facet_area[i] /=2.;
  }
  
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

  facet2vertex_average<Vector>( facet_norm, vertex_norm);
  for (int i=0;i<vertex_num;i++){
    vertex_norm[i] = vertex_norm[i] / CGAL::sqrt(vertex_norm[i] * vertex_norm[i]);
  }
}


void TriMesh::update_base(){//base update halfedge, facet, vertex.

  update_halfedge();

  update_facet();

  update_vertex();
};


void TriMesh::update_vertex_localchart(){
  vertex_LC[0] = Vec_Fun(vertex_num);
  vertex_LC[1] = Vec_Fun(vertex_num);
  
  for (int i=0;i<vertex_num;i++) localchart(vertex_LC[0][i],vertex_LC[1][i], vertex_norm[i]);
}

void TriMesh::update_facet_localchart(){
  facet_LC[0] = Vec_Fun(facet_num);
  facet_LC[1] = Vec_Fun(facet_num);

  for (int i=0;i<facet_num;i++) localchart(facet_LC[0][i],facet_LC[1][i], facet_norm[i]);
}


void TriMesh::update_facet_curvature(){

}
