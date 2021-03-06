/*
  FILE: TriMesh_Curv.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implements the curvature computation 
  of class TriMesh.
  
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

#include "meshtk/TriMesh.hh"

#include <fstream>

namespace meshtk{

  unsigned TriMesh::attribute_allocate(unsigned item, unsigned type, void * init_val){
    int n;

    if (item == MESHTK_VERTEX) n = vertex_num;
    else if (item == MESHTK_FACET) n = facet_num;
    else if (item == MESHTK_HALFEDGE) n = halfedge_num;
    else { std::cerr << "Attribute function: item indice is not correct, use MESHTK_VERTEX, MESHTK_FACET, MESHTK_HALFEDGE" << std::endl; exit(1);}

    if (type == MESHTK_SCALAR) {
      ScalarFunction *v;
      if (init_val !=NULL)
	v = new ScalarFunction(n, *((double *) init_val));
      else 
	v = new ScalarFunction(n);
      attribute[set_attribute_id] = v;
    }
    else if (type == MESHTK_VECTOR) {
      VectorFunction *v;
      if (init_val !=NULL)
	v = new VectorFunction(n, *((Vector *) init_val));
      else 
	v = new VectorFunction(n);
      attribute[set_attribute_id] = v;
    }
    else if (type == MESHTK_BOOLEAN) {
      BooleanFunction *v;
      if (init_val !=NULL)
	v = new BooleanFunction(n, *((bool *) init_val));
      else 
	v = new BooleanFunction(n);
      attribute[set_attribute_id] = v;
    }
    else { std::cerr << "Attribute function: type indice is not correct, use MESHTK_SCALAR, MESHTK_VECTOR, MESHTK_BOOLEAN" << std::endl; exit(1);}

    return set_attribute_id++;
  }

  void TriMesh::attribute_delete(unsigned indice, unsigned type){
    if (type == MESHTK_SCALAR) {
      ((ScalarFunction *) attribute[indice])->clear();
      delete (ScalarFunction *) attribute[indice];
    }
    else if (type == MESHTK_VECTOR) {
      ((VectorFunction *) attribute[indice])->clear();
      delete (VectorFunction *) attribute[indice];
    }
    else if (type == MESHTK_BOOLEAN) {
      ((BooleanFunction *) attribute[indice])->clear();
      delete (BooleanFunction *) attribute[indice];
    }
    else { std::cerr << "Warning: Attribute function type indice is not correct, use MESHTK_SCALAR, MESHTK_VECTOR, MESHTK_BOOLEAN" << std::endl; return; }    
    

    attribute.erase(indice);

  }

  void TriMesh::attribute_print(unsigned indice, unsigned type, std::string filename) {

    filename.append(".attr");
    std::ofstream file;
    file.open(filename.c_str());
    
    if (type == MESHTK_SCALAR) {
      ScalarFunction *v = (ScalarFunction *) attribute[indice];
      int n = v->size();
      for (int i = 0; i< n; ++i) 
	file << (*v)[i] << '\n';
    }
    else if (type == MESHTK_VECTOR) {
      VectorFunction *v = (VectorFunction *) attribute[indice];
      int n = v->size();
      for (int i = 0; i< n; ++i) 
	file << (*v)[i] << '\n';
    }
    else if (type == MESHTK_BOOLEAN) {
      BooleanFunction *v = (BooleanFunction *) attribute[indice];
      int n = v->size();
      for (int i = 0; i< n; ++i) 
	file << (*v)[i] << '\n';
    }
    else { std::cerr << "Attribute function: type indice is not correct, use MESHTK_SCALAR, MESHTK_VECTOR, MESHTK_BOOLEAN" << std::endl; exit(1);}
    
    file.close();
  }
  
  void * TriMesh::attribute_extract(unsigned indice){
    return attribute[indice];
  }

}
