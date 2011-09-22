#ifndef _TRIMESH_HH_
#define _TRIMESH_HH_
#include <vector>
#include "mesh_topo.hh"

namespace meshtk {

  typedef std::vector<Halfedge_handle> ISHalfedge_list;
  typedef std::vector<Vertex_handle> ISVertex_list;
  typedef std::vector<Facet_handle> ISFacet_list;

  typedef std::vector<Vector> Vec_Fun;
  typedef std::vector<double> Scalar_Fun;
  typedef std::vector<bool> Bool_Fun;// For marker.


  struct Curvature {
    Curvature (const double pc0, const double pc1, const double hc, const double kc) :PC0(pc0), PC1(pc1), hcurv(hc), kcurv(kc){}
    const double PC0, PC1, hcurv, kcurv;
    /*
    void operator= (const Curvature p){
      PC0 = p.PC0;
      PC1 = p.PC1;
      hcurv = p.hcurv;
      kcurv = p.hcurv;
    }
    */
  };

  class TriMesh {//topo ref system //index system of items of polyhedron

    Polyhedron P;

    double coord_min_x, coord_min_y, coord_min_z;
    double coord_max_x, coord_max_y, coord_max_z;


    Vec_Fun vertex_norm;//normal vector
    Vec_Fun facet_norm;


    //HalfedgeIS_list EI;
    //HalfedgeIS_list HI;

    //VertexIS_list VI;
    //FacetIS_list FI;

    ISHalfedge_list IH; // 
    ISVertex_list IV;
    ISFacet_list IF;

    //item attributes
    Vec_Fun halfedge_vec;


    Vec_Fun vertex_LC[2];//LC: local chart
    Vec_Fun facet_LC[2];
  
    Scalar_Fun vertex_CT[3];//CT: cuvature tensor [ 0 1 ]
    Scalar_Fun facet_CT[3];//                     [ 1 2 ]

    Scalar_Fun vertex_avg_len;
    
    Scalar_Fun vertex_area;
    Scalar_Fun facet_area;

    double avg_edge_len;


    double update_halfedge();
    double update_facet();
    void update_vertex();

    void update_facet_localchart();
    void update_vertex_localchart();


    void update_facet_curvature();
    void update_vertex_curvature();

  public:  

    double total_area;
    int halfedge_num;
    int vertex_num;
    int facet_num;


    Scalar_Fun vertex_PC0;//PC: principle cuvature
    Scalar_Fun vertex_PC1;
    Scalar_Fun facet_PC0;
    Scalar_Fun facet_PC1;

    Scalar_Fun vertex_hcurv;//mean curvature
    Scalar_Fun vertex_kcurv;//Gaussian curvature

    Scalar_Fun facet_hcurv;//mean curvature
    Scalar_Fun facet_kcurv;//Gaussian curvature


    void read(std::string, std::string);
    void write(std::string, std::string);
    void init_index();
    void update_base();
    void update_curvature();


    template <class T>
    void facet2vertex_area_average(std::vector<T> &f, std::vector<T> &v, T zero){
    };

    template <class T>
    void facet2vertex_point_average(std::vector<T> &f, std::vector<T> &v, T zero){
      double sigma = avg_edge_len;

      for (int i=0;i<vertex_num;i++){
	Vector tmp;
	double scale, total_scale=0;
	HV_circulator hv=IV[i]->vertex_begin();
	T vec = zero;

	sigma = vertex_avg_len[i] *2. /3.;
	
	do{
	  if (hv->facet()==NULL) continue;

	  tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
	  scale = CGAL::sqrt(tmp * tmp) /3.;		    
	  scale = std::exp(- scale * scale / (2 * sigma * sigma));	    

	  total_scale += scale;
	  vec = vec + scale * f[hv->facet()->index];

	} while (++ hv!=IV[i]->vertex_begin());


	v[i] = vec/total_scale;
      }
    };

    template <class T>
    void vertex2facet_average(std::vector<T> &v, std::vector<T> &f){
      for (int i=0;i<facet_num;i++){
	// ...
      }
    };






    friend class MeshPainter;
  

  };

}
#endif /* _TRIMESH_HH_ */
