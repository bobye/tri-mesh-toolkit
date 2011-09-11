#include <vector>
#include "mesh_topo.h"

void halfedgeUpdate(Polyhedron_IS & PI){
  
  int n=PI.P.size_of_halfedges();
  Halfedge_handle h;
  for (int i=0;i<n;i++){
    h=PI.IH[i];
    PI.halfedge_vec[i] = h->vertex()->point() - h->opposite()->vertex()->point();
  }
  
}


void facetUpdate(Polyhedron_IS & PI){
  int n=PI.P.size_of_facets();
  Vector normal;
  Halfedge_handle h;

  for (int i=0;i<n;i++){
    h=PI.IF[i]->halfedge();
    normal = CGAL::cross_product(halfedge_vec[PI.HI[h]], halfedge_vec[PI.HI[h->next()]]);
    PI.facet_area[i] = normal * normal ;
    PI.facet_norm[i] = normal / CGAL::sqrt(PI.facet_area[i]); 
    PI.facet_area[i] /=2.;
  }
}


void vertexUpdate(Polyhedron_IS & PI){
  
}


