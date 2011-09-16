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
  Vec_Fun facet_norm;
  Vec_Fun vertex_norm;

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
  void init_vertex_localchart(Vec_Fun *);
  void init_facet_localchart(Vec_Fun *);



};

#endif /* _TRIMESH_H_ */
