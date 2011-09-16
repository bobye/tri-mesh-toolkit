#ifndef _TRIMESH_H_
#define _TRIMESH_H_
#include <vector>
#include "mesh_topo.h"

typedef std::vector<Halfedge_handle> ISHalfedge_list;
typedef std::vector<Vertex_handle> ISVertex_list;
typedef std::vector<Facet_handle> ISFacet_list;

typedef std::vector<Vector> Vec_Fun;
typedef std::vector<double> Scalar_Fun;


class TriMesh {//topo ref system //index system of items of polyhedron
  Polyhedron P;

  //HalfedgeIS_list EI;
  //HalfedgeIS_list HI;

  //VertexIS_list VI;
  //FacetIS_list FI;

  ISHalfedge_list IH; // 
  ISVertex_list IV;
  ISFacet_list IF;

  //item attributes
  Vec_Fun halfedge_vec;

  Vec_Fun vertex_norm;
  Vec_Fun facet_norm;

  Vec_Fun vertex_LC[2];
  Vec_Fun facet_LC[2];


  Scalar_Fun facet_area;

  double avg_edge_len;


  double update_halfedge();
  double update_facet();
  void update_vertex();

public:  
  TriMesh (){};

  double total_area;
  int edge_num;
  int vertex_num;
  int facet_num;


  void read(std::string, std::string);
  void write(std::string, std::string);
  void init_index();
  void update_base();

  template <class T>
  void facet2vertex_average(std::vector<T> &f, std::vector<T> &v){
    double sigma = 2 * avg_edge_len;

    for (int i=0;i<vertex_num;i++){
      Vector tmp;
      double scale, total_scale=0;
      HV_circulator hv=IV[i]->vertex_begin();
      T vec;
      
      do {
	tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
	scale = CGAL::sqrt(tmp * tmp);
	total_scale += scale = std::exp(- scale * scale / (sigma * sigma));
	vec = vec + scale * f[hv->facet()->index];	
      }while (++hv!=IV[i]->vertex_begin());
      v[i] = vec/total_scale;
    }
  };

  template <class T>
  void vertex2facet_average(std::vector<T> &v, std::vector<T> &f){
    for (int i=0;i<facet_num;i++){
      // ...
    }
  };

  void update_vertex_localchart();
  void update_facet_localchart();
  void update_facet_curvature();




};

#endif /* _TRIMESH_H_ */
