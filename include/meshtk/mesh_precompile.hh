/*
  FILE: meshtk/MeshPrecompile.hh This file is part of tri-mesh-toolkit.
  It is a C++ header file which defines global constant for precompile
  
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
#ifndef _MESH_PRECOMPILE_HH_
#define _MESH_PRECOMPILE_HH_

#define MESHTK_VERSION                 ("0.4.1")



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
#define MESHTK_VERTEX_NORM              32

#define MESHTK_FACET_NORM               33
// Boolean section: 64 - 95             
#define MESHTK_VERTEX_SALIENT           64
#define MESHTK_VERTEX_SALIENT_SUP       65
#define MESHTK_VERTEX_SALIENT_INF       66

// Point section: 96 - 127
#define MESHTK_VERTEX_COORD             96






#define MESHTK_PI                       (3.1415926535898)
#define MESHTK_SCALAR_TYPE              double

#endif /* _MESH_PRECOMPILE_HH_ */
