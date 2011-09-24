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

    // just the coordinate buffer for vertex and facet
    PointFunction vertex_coord;
    PointFunction facet_coord;
    
    // Salient vertex mark
    BooleanFunction vertex_salient;
    BooleanFunction vertex_salient_sup;
    BooleanFunction vertex_salient_inf;
    /////////////////////////////////////////////////////////////////////
    // private routines to update item attributes starts here

  public:

    /////////////////////////////////////////////////////////////////////////////
    // public functions used in top interface.

    // build connection with lower CGAL layer, should be called
    // immediatelly after loading the mesh, update: IH, IV, IF
    // and [vertex, facet, halfedge]_handle->index 
    void init_index();    

    // push vertex_coord to CGAL Polyhedron data mesh
    void restore_coord();

    // Apply Gaussian smooth operator over entire mesh
    // The new point is an averaging of one ring configuration with unnormalized
    //    weight = exp ( - distance ^2 / sigma ^2 )
    // where distance is the attached edge length, and the argument coeff is given by
    //    sigma = coeff * avg_edge_len
    void gaussian_smooth(double );

    // Update mark for salient vertex, the argument are taken as num of smooth 
    // iteration (>0 preconditioned)
    void update_vertex_salient(int , int pre_inter = 1);




  };
}


#endif /* _DYNAMICTRIMES_HH_ */
