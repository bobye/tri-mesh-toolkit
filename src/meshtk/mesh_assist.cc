/*
  FILE: mesh_assist.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implements some assistant functions.
  
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


#include "meshtk/mesh_topo.hh"
#include <ctime>


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


double Heron_formula(double l1, double l2, double l3) {
  double p = (l1+l2+l3)/2;
  return std::sqrt((p-l1) * (p-l2) * (p-l3) * p);
}

clock_t start, end;
void clock_start(std::string description) {
  start = clock();
  std::cout << description <<" ... " << std::flush; 
}

void clock_end() {
  end = clock();
  std::cout << "\t[done] " << (static_cast<double> (end) - static_cast<double> (start)) / CLOCKS_PER_SEC <<" seconds" << std::endl;
}



