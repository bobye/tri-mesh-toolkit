#include <vector>
#include "mesh_topo.h"

double halfedgeUpdate(Polyhedron_IS & PI){
  
  int n=PI.edge_num;
  double avg_len=0;
  Halfedge_handle h;
  for (int i=0;i<2*n;i++) 
    if ((h=PI.IH[i])!=NULL){    
      PI.halfedge_vec[i] = h->vertex()->point() - h->opposite()->vertex()->point();
      avg_len += CGAL::sqrt(PI.halfedge_vec[i] * PI.halfedge_vec[i]);
    }
  return avg_len/(2*n-PI.P.size_of_border_edges());
}


double facetUpdate(Polyhedron_IS & PI){
  int n=PI.P.size_of_facets();
  double total_area=0;
  Vector normal;
  Halfedge_handle h;

  for (int i=0;i<n;i++){
    h=PI.IF[i]->halfedge();
    normal = CGAL::cross_product(PI.halfedge_vec[h->index], PI.halfedge_vec[h->next()->index]);
    total_area += PI.facet_area[i] = normal * normal ;
    PI.facet_norm[i] = normal / CGAL::sqrt(PI.facet_area[i]); 
    PI.facet_area[i] /=2.;
  }
  return total_area;
}


void vertexUpdate(Polyhedron_IS & PI){
  int n=PI.P.size_of_vertices();
  double sigma = 2 * PI.avg_edge_len;

  for (int i=0;i<n;i++){
    Vector normal(0,0,0), tmp;
    double scale;
    HV_circulator hv=PI.IV[i]->vertex_begin();
    do{
      tmp = (-PI.halfedge_vec[hv->index]+PI.halfedge_vec[hv->next()->index]);
      scale = CGAL::sqrt(tmp * tmp);
      normal = normal + std::exp(- scale * scale / (sigma *sigma)) * PI.facet_norm[hv->facet()->index];
    }while (++hv!=PI.IV[i]->vertex_begin());
    PI.vertex_norm[i] = normal / CGAL::sqrt(normal * normal);
  }
}

