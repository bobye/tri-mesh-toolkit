/*
  FILE: meshtk/KeyPoint.hh This file is part of MeshTK.
  It is a C++ header file which defines the interface of class KeyPoint. 
  
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
#ifndef _KEYPOINT_HH_
#define _KEYPOINT_HH_

#include <map>
#include "MeshAttributeType.hh"

#define MESHTK_SIFT_BINS_NUMBER (36)

namespace meshtk {


  struct KeyPoint {

  public:  

    int index; //reference in mesh vertex
    double scale; // scale distance in mesh surface
    double magnitude; // the extreme of DoH

    ScalarNeighborFunction vertex_neighbor;
    NeighborIndex facet_neighbor;

    double histogram[MESHTK_SIFT_BINS_NUMBER];// descriptor, 36 bins [10*i-5, 10*i+5), i=0,...,35


    KeyPoint(int i, double s, double m) :
      index(i), scale(s), magnitude(m) {
    }
    
  
  };


  struct SIFTPoint {
  };



}
#endif /* _KEYPOINT_HH_ */
