#ifndef _PERELTOPO_H
#define _PERELTOPO_H

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
  CGAL::Color color;//
};


template <class Refs>
class Perel_face : public CGAL::HalfedgeDS_face_base<Refs> {
public:
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


typedef std::map<Halfedge_handle, int, Handle_comparison<Halfedge_handle> > HalfedgeIS_list;
typedef std::map<Vertex_handle, int, Handle_comparison<Vertex_handle> >     VertexIS_list;
typedef std::map<Facet_handle, int, Handle_comparison<Facet_handle> >       FacetIS_list;
typedef std::vector<Halfedge_handle> ISHalfedge_list;
typedef std::vector<Vertex_handle> ISVertex_list;
typedef std::vector<Facet_handle> ISFacet_list;

struct Polyhedron_IS {//topo ref system
  Polyhedron P;

  HalfedgeIS_list EI;
  HalfedgeIS_list HI;

  VertexIS_list VI;
  FacetIS_list FI;

  ISHalfedge_list IH; // 
  ISVertex_list IV;
  ISFacet_list IF;
  
  
};

struct perelGuassianInit {
  template <class Vertex>
  void operator()(Vertex &v){
  }
};

struct perelLaplaceKInit {
  template <class Vertex>
  void operator()(Vertex &v){
  }
};



struct Polyhedron_Init{
    template <class Polyhedron_IS>    
    void operator()(Polyhedron_IS& PI){
      //        std::for_each( P.vertices_begin(), P.vertices_end(), perelGuassianInit());
      //  std::for_each( P.vertices_begin(), P.vertices_end(), perelLaplaceKInit());
      
      int n=(PI.P.size_of_halfedges()+PI.P.size_of_border_edges())/2;
      PI.IH = ISHalfedge_list(2*n);
      PI.IV =  ISVertex_list(PI.P.size_of_vertices());
      PI.IF =  ISFacet_list(PI.P.size_of_facets());

      int index_count=0;
      for(Vertex_iterator vitr= PI.P.vertices_begin();vitr!= PI.P.vertices_end();
	  PI.IV[index_count]=vitr, PI.VI[vitr] = index_count++, vitr++);
      index_count=0;
      for(Edge_iterator eitr= PI.P.edges_begin();eitr!= PI.P.edges_end();
	  PI.IH[index_count]=eitr, PI.IH[index_count+n]=eitr->opposite(),
	    PI.EI[eitr] = PI.EI[eitr->opposite()] = index_count, 
	    PI.HI[eitr] = index_count, PI.HI[eitr->opposite()] = index_count +n,
	    index_count++, eitr++);
      index_count=0;
      for(Facet_iterator fitr= PI.P.facets_begin(); fitr!= PI.P.facets_end(); PI.IF[index_count]=fitr, PI.FI[fitr]=index_count++,fitr++);      
    };
};



#endif
