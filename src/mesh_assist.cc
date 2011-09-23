/*
  FILE: mesh_assist.cc This file is part of MeshTK.
  It is a C++ source file which implements some assistant functions.
  
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


#include "meshtk/mesh_topo.hh"


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








