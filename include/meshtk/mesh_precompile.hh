/*
  FILE: meshtk/MeshPrecompile.hh This file is part of MeshTK.
  It is a C++ header file which defines global constant for precompile
  
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
#ifndef _MESH_PRECOMPILE_HH_
#define _MESH_PRECOMPILE_HH_

#define MESHTK_USER_ATTRIBUTE_START (128)
///////////////////////////////////////////////////////////////
//Following are used in function TriMesh::attribute_allocate()
// 1st argument
#define MESHTK_VERTEX                   0
#define MESHTK_FACET                    1
#define MESHTK_HALFEDGE                 2
// 2nd argument
#define MESHTK_SCALAR                   3
#define MESHTK_VECTOR                   4
#define MESHTK_BOOLEAN                  5
#define MESHTK_POINT                    6
/////////////////////////////////////////////////////////////////
// Following are registration number of inherent mesh attributes
// Scalar section: 0-31
#define MESHTK_VERTEX_PC0               0
#define MESHTK_VERTEX_PC1               1
#define MESHTK_VERTEX_HCURV             2
#define MESHTK_VERTEX_KCURV             3

#define MESHTK_FACET_PC0                4
#define MESHTK_FACET_PC1                5
#define MESHTK_FACET_HCURV              6
#define MESHTK_FACET_KCURV              7
// Vector section: 32-63


#endif /* _MESH_PRECOMPILE_HH_ */
