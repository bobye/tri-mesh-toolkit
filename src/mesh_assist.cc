#include "meshtk/mesh_topo.hh"
#include "meshtk/TriMesh.hh"


void localchart(Vector& u, Vector& v, Vector norm){
  Vector tmp;
  norm = norm / CGAL::sqrt( norm * norm);

  do {
    tmp = Vector (rand(), rand(), rand());
    u = tmp - (tmp * norm) * norm; 
  }while ( u*u < 1e-2 );

  u = u/ CGAL::sqrt(u * u);
  v = CGAL::cross_product(norm, u);
  v = v/ CGAL::sqrt(v * v);

}



void localcoord(Vector proj, Vector u, Vector v, double coord[2]){
  coord[0]=proj * u;
  coord[1]=proj * v;
}








