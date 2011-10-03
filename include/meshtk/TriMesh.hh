/*
  FILE: meshtk/TriMesh.hh This file is part of MeshTK.
  It is a C++ header file which defines the interface of class TriMesh. 
  
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

#ifndef _TRIMESH_HH_
#define _TRIMESH_HH_
#include <vector>
#include <set>
#include "mesh_topo.hh"
#include "mesh_precompile.hh"

namespace meshtk {

  // ISHalfedge_list represents Index System of Halfedges by std::vector. 
  // After initialization of class TriMesh by
  //    ISHalfedge_list IH;
  // One could use data member IH to reference a halfedge_handle instance
  // For example:
  //    Vertex_handle v = IH[i]->vertex(); //refer the vertex attached with
  //                                       //the i-th halfedge
  // ISVertex_list and ISFacet_list can be used in the same way.
  typedef std::vector<Halfedge_handle> ISHalfedgeList;
  typedef std::vector<Vertex_handle> ISVertexList;
  typedef std::vector<Facet_handle> ISFacetList;

  // Feature type defined from mesh, VectorFunction represents 3D vector function
  // define over mesh domain, with respect to halfedges, vertices or facets.
  // ScalarFunction and BooleanFunction correspond to scalar and boolean
  // function defined over mesh domain.
  typedef std::vector<Vector> VectorFunction; // 
  typedef std::vector<double> ScalarFunction; // displayed by color ramper
  typedef std::vector<double> BooleanFunction;// displayed by point marker.

  // Curvature data type
  struct Curvature {
  public:
    Curvature (const double pc0, const double pc1, const double hc, const double kc) 
 :principle_curv0(pc0), principle_curv1(pc1), mean_curv(hc), gaussian_curv(kc){}
    const double principle_curv0, principle_curv1, mean_curv, gaussian_curv;
    
    static Curvature tensor_compute(double, double, double);
  };


  // class TriMesh is the main part of MeshTK, which manages almost all mesh processing functions
  // This class may not provide access to individual vertex or facet of mesh. User only manipulate
  // mesh globally, like compute normal, estimate curvature, denoise or smooth mesh, etc. and 
  // manage data by returning a pointer of array.
  class TriMesh {/* One may use TriMesh for static mesh analysis, if the mesh is dynamic, such 
		    as smoothing, denoising, deformation and remeshing,  please consider the 
		    inheritant class called DynamicTriMesh  */
		 
  protected:
    // CGAL Polyhedron data
    Polyhedron P;

    // the bounding box of mesh
    double coordinate_min_x, coordinate_min_y, coordinate_min_z;
    double coordinate_max_x, coordinate_max_y, coordinate_max_z;

    // normal function defined on mesh domain w.r.t vertex and facet
    VectorFunction vertex_norm;
    VectorFunction facet_norm;

    // index system provide random access to object handles
    ISHalfedgeList IH; // 
    ISVertexList IV;
    ISFacetList IF;

    /////////////////////////////////////////////////////////////////////
    // inherent item attributes starts here
    
    // halfedge vector
    VectorFunction halfedge_vec;    

    // local chart(coordinate axis) used in tangent projection 
    // to achieve local parametrization 
    VectorFunction vertex_LC[2];//LC: local chart
    VectorFunction facet_LC[2];
  
    // locate curvature tensor with respect to *_LC (local chart)
    // 2x2 symmetric matrix 
    ScalarFunction vertex_CT[3];//CT: cuvature tensor [ 0 1 ]
    ScalarFunction facet_CT[3];//                     [ 1 2 ]

    // averge edge length attached to per-vertex
    ScalarFunction vertex_avg_len;
    
    // area dominated, vertex dominant area is simply compute by averaging
    // surrounding facets.
    ScalarFunction vertex_area;
    ScalarFunction facet_area;

    // two principle curvature, while PC0 is the larger one
    ScalarFunction vertex_PC0;//PC: principle cuvature
    ScalarFunction vertex_PC1;
    ScalarFunction facet_PC0;
    ScalarFunction facet_PC1;

    // hcurv = ( PC0 + PC1 )/2, kcurv = PC0 * PC1
    ScalarFunction vertex_hcurv;//mean curvature
    ScalarFunction vertex_kcurv;//Gaussian curvature
    ScalarFunction facet_hcurv;
    ScalarFunction facet_kcurv;


    // neighbor vertices indices, which store the a neighborhood of vertex
    // the indices of neighbor are stored in set with an iterator access
    std::vector<std::set<int> > vertex_neighbor;
    //std::vector<std::set<int> > facet_neighbor;

    double avg_edge_len;//average edge length globally


    /////////////////////////////////////////////////////////////////////
    // private routines to update item attributes starts here

    double update_halfedge();// to update: halfedge_vec, avg_edge_len    
    double update_facet();// to update: facet_norm, facet_area
    void update_vertex();// to update: vertex_norm, vertex_area, vertex_avg_len

    void update_facet_localchart();// to update local chart setting of facet - facet_LC[]
    void update_vertex_localchart();// to update local chart setting of vertex - vertex_LC[]

    void update_facet_curvature();// to update: facet_CT[], facet_PC*, facet_?curv
    void update_vertex_curvature();// to update: vertex_CT[], vertex_PC*, vertex_?curv


    ///////////////////////////////////////////////////////////////////////
    // Map register number to reference of functions define over mesh domain. 
    // Register number is unsigned integer and have some predefined ones for inherent item
    // attributes like coordinate, normal and curvature. 
    // Users may request additional register numbers for use of creating/accessing the 
    // attribute functions
    std::map<unsigned, void *> attribute; 
    // available register number to set
    // scalar inherent attribute range from 0-31 and vector inherent attribute range from 32-63
    // 64-127 is left for extension. user customed number started at 128.
    unsigned set_attribute_id; 

  public:  
    TriMesh();
    ~TriMesh();
    ////////////////////////////////////////////////////////////////////////////
    double total_area; // total area of mesh surface

    // The mesh structures in general are not permitted to changed.
    int halfedge_num;
    int vertex_num;
    int facet_num;


    
    /////////////////////////////////////////////////////////////////////////////
    // public functions used in top interface.

    // read mesh with format specification, for example 
    //   M.read("input", "off"); 
    // would load input.off into the TriMesh instance M
    void read(std::string , std::string);
    // like read( , ), just output mesh with specific format
    void write(std::string, std::string);

    // build connection with lower CGAL layer, should be called
    // immediatelly after loading the mesh, update: IH, IV, IF
    // and [vertex, facet, halfedge]_handle->index 
    virtual void init_index();

    // update base attributes of mesh, namely three private routines:
    //  update_halfedge(), update_facet(), update_vertex(); 
    void update_base();

    // update curvature attributes of mesh, update: 
    // update_[facet, vertex]_localchart() and update_[facet, vertex]_curvature()
    void update_curvature();

    /**************************************************************************/

    // update vertices neighbor, argument coeff is to take all *-ring within a Euclidean distance
    // return the average number of neighbor vertices associated
    double update_vertex_neighbor(double );//to update: vertex_neighbor[][]

    // guassian smooth an attribute function 
    // input: v0
    // output: v1
    template <class T>
    void gaussian_smooth_vertex(double coeff, std::vector<T> &v0, std::vector<T> &v1, T zero){
      // preconditioned with vertex_neighbor 
      // update_vertex_neighbor(3 * coeff);
      double sigma = coeff * avg_edge_len;


      for (int i = 0; i < vertex_num; i++){
	T vec = zero;
	Vector tmp;
	double scale, total_scale = 0;
      

	for (std::set<int>::iterator it = vertex_neighbor[i].begin();
	     it != vertex_neighbor[i].end(); ++it) {
	  tmp = IV[*it]->point() - IV[i]->point();
	  scale = CGAL::sqrt(tmp * tmp);
	  scale = vertex_area[*it] * std::exp( - (scale * scale) / (2 * sigma * sigma));
	  total_scale += scale;
	  vec = vec + scale * v0[*it];
	}

	v1[i] = vec / total_scale;
      
      }
      
    };

    // smooth with prescribed coefficient, normalized precondition
    template <class T>
    void smooth_vertex(std::vector<std::map<int, double> > &coeff, 
				std::vector<T> &v0, std::vector<T> &v1, T zero){
      for (int i = 0; i < vertex_num; i++){
	v1[i] = zero;	

	for (std::map<int, double>::iterator it = coeff[i].begin();
	     it != coeff[i].end(); ++it) {
	  v1[i] = v1[i] + (*it).second * v0[(*it).first];
	}
      
      }

    }
    
    // The following procedure is SIFT keypoint detection for scalar 
    // function on static manifold mesh domain. The input is scalar
    // function, the keypoints detected are given by boolean function
    int detect_vertex_keypoint(ScalarFunction &, BooleanFunction &, int, int pre_iter=1);
    

    /**************************************************************************/
    // allocate memory for attribute function
    unsigned attribute_allocate(unsigned, unsigned);
    // return reference of attribute function by register number
    void *attribute_extract(unsigned );
    void attribute_delete(unsigned, unsigned);


    /**************************************************************************/
    
    // Below are weighted average template functions, which will 
    // convert facet attribute to vertex attribute or vice versa
    template <class T>
    void facet2vertex_area_average(std::vector<T> &f, std::vector<T> &v, T zero){
    };// not implemented yet

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
	  scale = facet_area[hv->facet()->index] * std::exp(- scale * scale / (2 * sigma * sigma));	    

	  total_scale += scale;
	  vec = vec + scale * f[hv->facet()->index];

	} while (++ hv!=IV[i]->vertex_begin());


	v[i] = vec/total_scale;
      }
    };

    template <class T>
    void vertex2facet_average(std::vector<T> &v, std::vector<T> &f){
      for (int i=0;i<facet_num;i++){
	// ... still not implemented 
      }
    };



    template <class T>
    void vertex2vertex_average(std::vector<T> &s, std::vector<T> &t){
      // not implemented yet
    };


    // some protected variable are accessible in MeshPainter and its derivative
    friend class MeshPainter;
  

  };

}
#endif /* _TRIMESH_HH_ */
