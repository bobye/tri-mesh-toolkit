/*
  FILE: meshtk/mesh_assist.hh This file is part of tri-mesh-toolkit.
  It is a C++ header file which declares some assistant functions. 
  
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

#ifndef _MESH_ASSIST_HH_
#define _MESH_ASSIST_HH_

#include "mesh_topo.hh"
#include <string>

extern void localchart(Vector &, Vector &, Vector);
extern void localcoord(Vector, Vector, Vector, double*);

extern double Heron_formula(double, double, double);

template <class T> 
struct index_cmp {
  index_cmp(const T arr): arr(arr) {}
  bool operator() (const size_t a, const size_t b) const {
    return arr[a] < arr[b];
  }
  const T arr;
};

extern void clock_start(std::string );
extern void clock_end();
#endif /* _MESH_ASSIST_HH_ */
