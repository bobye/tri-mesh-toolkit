#include "mesh_topo.h"

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


double prin_curv(double e, double f, double g, double &k1, double &k2){
  double H= (e+g)/2;
  double G=std::sqrt(4*f*f-(e-g)*(e-g))/2;
  k1=H+G;
  k2=H-G;
  return H;
}




