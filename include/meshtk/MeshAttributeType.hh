/*
  FILE: meshtk/MeshAttributeType.hh This file is part of tri-mesh-toolkit.
  It is a C++ header file which declares mesh attribute types
  
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



#ifndef _MESHATTRIBUTETYPE_HH_
#define _MESHATTRIBUTETYPE_HH_

#include "mesh_precompile.hh"

namespace meshtk {

  

  // ISHalfedge_list represents Index System of Halfedges by std::vector. 
  // After initialization of class TriMesh by
  //    ISHalfedge_list IH;
  // One could use data member IH to reference a halfedge_handle instance
  // For example:
  //    Vertex_handle v = IH[i]->vertex(); //refer the vertex attached with
  //                                       //the i-th halfedge
  // ISVertex_list and ISFacet_list can be used in the same way.
  typedef std::vector<Halfedge_handle> ISHalfedgeList;
  typedef std::vector<Vertex_handle> ISVertexList;
  typedef std::vector<Facet_handle> ISFacetList;

  // Feature type defined from mesh, VectorFunction represents 3D vector function
  // define over mesh domain, with respect to halfedges, vertices or facets.
  // ScalarFunction and BooleanFunction correspond to scalar and boolean
  // function defined over mesh domain.
  typedef std::vector<Vector> VectorFunction; // 
  typedef std::vector<MESHTK_SCALAR_TYPE> ScalarFunction; // displayed by color ramper
  typedef std::vector<bool> BooleanFunction;// displayed by point marker.

  // type used to load and restore coordinate of points
  typedef std::vector<Point> PointFunction;
  typedef std::vector<PointFunction> PointFuncSequence;

  // neighbor faces and vertices
  typedef std::map<int, double> ScalarNeighborFunction;
  typedef std::set<int> NeighborIndex;

  // Curvature data type
  struct Curvature {
  public:
    Curvature (const double pc0, const double pc1, const double hc, const double kc) 
 :principle_curv0(pc0), principle_curv1(pc1), mean_curv(hc), gaussian_curv(kc){}
    const double principle_curv0, principle_curv1, mean_curv, gaussian_curv;
    
    static Curvature tensor_compute(double, double, double);
  };

}


#endif /* _MESHATTRIBUTETYPE_HH_ */
