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
  typedef std::vector<PointFunction>           PointFuncSequence;
  typedef std::vector<VectorFunction>          VectorFuncSequence;
  typedef std::vector<int>                     LabelFunction;

  class ManifoldTriMesh : public TriMesh {
    PointFuncSequence   vertex_coord_seq;
    VectorFuncSequence  vertex_rotate_seq;
    VectorFuncSequence  vertex_rotate_exp;

    ScalarFunction halfedge_cot;    

  public:
    /////////////////////////////////////////////////////////////////////////////
    // public functions used in top interface.
    LabelFunction vertex_label;
    std::vector<float> label_color;

    void update_base();
    void load_sequence(std::string );
    void load_examples(std::string );

    void PETSc_assemble_graphcut();
    void PETSc_export_graphcut_vectors(std::string );
    
    void compute_rotate_sequence();
    void print_rotate_sequence(std::string);
    
    void load_proxy_bone(std::string);
  };
}
#endif /* _MANIFOLDTRIMESH_HH_ */

