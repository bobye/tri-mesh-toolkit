#ifndef _MESH_TOPO_H
#define _MESH_TOPO_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>


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

template<class Handle>
struct Handle_comparison
{
  bool operator() (const Handle& h1, const Handle& h2) const
  {
    return &(*h1) < &(*h2);
  }
}; 

//map is not efficient
//typedef std::map<Halfedge_handle, int, Handle_comparison<Halfedge_handle> > HalfedgeIS_list;
//typedef std::map<Vertex_handle, int, Handle_comparison<Vertex_handle> >     VertexIS_list;
//typedef std::map<Facet_handle, int, Handle_comparison<Facet_handle> >       FacetIS_list;

//random access to items
typedef std::vector<Halfedge_handle> ISHalfedge_list;
typedef std::vector<Vertex_handle> ISVertex_list;
typedef std::vector<Facet_handle> ISFacet_list;

struct Polyhedron_IS {//topo ref system //index system of items of polyhedron
  Polyhedron P;

  //HalfedgeIS_list EI;
  //HalfedgeIS_list HI;

  //VertexIS_list VI;
  //FacetIS_list FI;

  ISHalfedge_list IH; // 
  ISVertex_list IV;
  ISFacet_list IF;

  //item attributes
  std::vector<Vector> halfedge_vec;
  std::vector<Vector> facet_norm;
  std::vector<Vector> vertex_norm;

  std::vector<double> facet_area;
  
};



struct Polyhedron_Init{
  template <class Polyhedron_IS>    
  void operator()(Polyhedron_IS& PI){
      
    int n=(PI.P.size_of_halfedges()+PI.P.size_of_border_edges())/2;
    PI.IH = ISHalfedge_list(2*n);
    PI.IV =  ISVertex_list(PI.P.size_of_vertices());
    PI.IF =  ISFacet_list(PI.P.size_of_facets());
      
    PI.halfedge_vec = std::vector<Vector> (2*n);
    PI.facet_norm   = std::vector<Vector> (PI.P.size_of_facets());
    PI.vertex_norm  = std::vector<Vector> (PI.P.size_of_vertices());
    PI.facet_area   = std::vector<double> (PI.P.size_of_facets());

    int index_count=0;
    for(Vertex_iterator vitr= PI.P.vertices_begin();vitr!= PI.P.vertices_end();
	PI.IV[index_count]=vitr, vitr->index = index_count++, vitr++);
    index_count=0;
    for(Edge_iterator eitr= PI.P.edges_begin();eitr!= PI.P.edges_end();
	PI.IH[index_count]=eitr, PI.IH[index_count+n]=eitr->opposite(),
	  //	  PI.EI[eitr] = PI.EI[eitr->opposite()] = index_count, 
	  eitr->index = index_count, eitr->opposite()->index = index_count +n,
	  index_count++, eitr++);
    index_count=0;
    for(Facet_iterator fitr= PI.P.facets_begin(); fitr!= PI.P.facets_end(); PI.IF[index_count]=fitr, fitr->index=index_count++,fitr++);      
  };

};



#endif
