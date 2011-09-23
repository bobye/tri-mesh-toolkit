/*
  FILE: meshtk/DynamicTriMesh.hh This file is part of MeshTK.
  It is a C++ header file which defines the interface of class DynamicTriMesh. 
  
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

#ifndef _DYNAMICTRIMES_HH_
#define _DYNAMICTRIMES_HH_

#include "TriMesh.hh"

namespace meshtk {
  // type used to load and restore coordinate of points
  typedef std::vector<Point> PointFunction;

  // DynamicTriMesh is inherited from TriMesh. It has ability to be applied some 
  // dynamic operations, such like smoothing, denoising, and remeshing.
  class DynamicTriMesh : public TriMesh {
    
    PointFunction vertex_coord;
    PointFunction facet_coord;
    
    /////////////////////////////////////////////////////////////////////
    // private routines to update item attributes starts here

    virtual double update_halfedge();// to update: halfedge_vec, avg_edge_len    
    virtual double update_facet();// to update: facet_norm, facet_area
    virtual void update_vertex();// to update: vertex_norm, vertex_area, vertex_avg_len

  public:

    /////////////////////////////////////////////////////////////////////////////
    // public functions used in top interface.
    
    // update base attributes of mesh, namely three private routines:
    //  update_halfedge(), update_facet(), update_vertex(); 
    virtual void update_base();
    
  }
}


#endif /* _DYNAMICTRIMES_HH_ */
