/*
  FILE: meshtk/mesh_assist.hh This file is part of MeshTK.
  It is a C++ header file which declares some assistant functions. 
  
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

#ifndef _MESH_ASSIST_HH_
#define _MESH_ASSIST_HH_

#include "mesh_topo.hh"

extern void localchart(Vector &, Vector &, Vector);
extern void localcoord(Vector, Vector, Vector, double*);


template <class T> 
struct index_cmp {
  index_cmp(const T arr): arr(arr) {}
  bool operator() (const size_t a, const size_t b) const {
    return arr[a] < arr[b];
  }
  const T arr;
};


#endif /* _MESH_ASSIST_HH_ */
