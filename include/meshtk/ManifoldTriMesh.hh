/*
  FILE: meshtk/ManifoldTriMesh.hh This file is part of tri-mesh-toolkit.
  It is a C++ header file which defines the interface of class ManifoldTriMesh. 
  
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
#ifndef _MANIFOLDTRIMESH_HH_
#define _MANIFOLDTRIMESH_HH_

#include "TriMesh.hh"
namespace meshtk {
  class ManifoldTriMesh : public TriMesh {
    PointFuncSequence vertex_coord_seq;
    
  public:
    /////////////////////////////////////////////////////////////////////////////
    // public functions used in top interface.

    // build connection with lower CGAL layer, should be called
    // immediatelly after loading the mesh, update: IH, IV, IF
    // and [vertex, facet, halfedge]_handle->index 
    void init_index();    
    
  };
}
#endif /* _MANIFOLDTRIMESH_HH_ */

