/*
  FILE: meshtk/mesh_topo.hh This file is part of MeshTK.
  It is a C++ header file which provides a wrapper for CGAL Polyhedron 
  data structure. 
  
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

#ifndef _MESH_TOPO_HH_
#define _MESH_TOPO_HH_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>



typedef CGAL::Simple_cartesian<double>                       Kernel;
//typedef CGAL::Cartesian<double>                              Kernel;
typedef Kernel::Vector_3                                     Vector;
//typedef Kerenl::Vector_2                                     Vector2D;
typedef Kernel::Point_3                                      Point;
//typedef Kernel::Point_2                                      Point2D;

template <class Refs, class T, class P>
class Perel_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
public:
  Perel_vertex() {} // repeat mandatory constructors
  Perel_vertex( const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {}
  //CGAL::Color color;//
  int index;
};


template <class Refs>
class Perel_face : public CGAL::HalfedgeDS_face_base<Refs> {
public:
  int index;
};

template <class Refs>
class Perel_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs> {
public:
  int index;
};

struct Perel_items : public CGAL::Polyhedron_items_3 {
  template <class Refs, class Traits>
  struct Vertex_wrapper {
    typedef typename Traits::Point_3  Point;
    typedef Perel_vertex<Refs, CGAL::Tag_true, Point> Vertex;
  };

  template <class Refs, class Traits>
  struct Face_wrapper {
    typedef Perel_face<Refs> Face;
  };

  template <class Refs, class Traits>
  struct Halfedge_wrapper {
    typedef Perel_halfedge<Refs> Halfedge;
  };

};






typedef CGAL::Polyhedron_3<Kernel, Perel_items>                 Polyhedron;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_handle                            Vertex_handle;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;

typedef Polyhedron::Halfedge                                 Halfedge;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Halfedge_iterator                        Halfedge_iterator;

typedef Polyhedron::Facet                                    Facet;
typedef Polyhedron::Facet_handle                             Facet_handle;
typedef Polyhedron::Facet_iterator                           Facet_iterator;

typedef Polyhedron::Edge_iterator                            Edge_iterator;


typedef Polyhedron::Halfedge_around_vertex_circulator        HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;






#endif
