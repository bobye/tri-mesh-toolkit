#ifndef _MESH_TOPO_H
#define _MESH_TOPO_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>



typedef CGAL::Simple_cartesian<double>                       Kernel;
//typedef CGAL::Cartesian<double>                              Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;


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


/* template<class Handle> */
/* struct Handle_comparison */
/* { */
/*   bool operator() (const Handle& h1, const Handle& h2) const */
/*   { */
/*     return &(*h1) < &(*h2); */
/*   } */
/* };  */

//map is not efficient
//typedef std::map<Halfedge_handle, int, Handle_comparison<Halfedge_handle> > HalfedgeIS_list;
//typedef std::map<Vertex_handle, int, Handle_comparison<Vertex_handle> >     VertexIS_list;
//typedef std::map<Facet_handle, int, Handle_comparison<Facet_handle> >       FacetIS_list;

//random access to items





#endif
