/*
  FILE: ManifoldTriMesh.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implemented class ManifoldTriMesh.
  
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

#include <iostream>
#include <fstream>
#include "meshtk/ManifoldTriMesh.hh"
#include "fsvd.hh"

namespace meshtk {  
  void ManifoldTriMesh::update_base() {
    TriMesh::update_base();
    halfedge_cot.resize(halfedge_num);
    Halfedge_handle h;
    for (int i=0; i<halfedge_num; ++i) 
      if ((h = IH[i])->facet() != NULL){
	Vector a = halfedge_vec[h->next()->index],
	  b = halfedge_vec[h->prev()->index];
	halfedge_cot[i] = a*b/ CGAL::sqrt(CGAL::cross_product(a,b).squared_length());
      } else {
	halfedge_cot[i] = 0;
      }
  }


  void ManifoldTriMesh::load_sequence(std::string filename) {
    filename.append(".pts");
    std::fstream fid; fid.open(filename.c_str());
    int n, m=vertex_coord_seq.size(); fid >> n;

    vertex_coord_seq.resize(m+n);
    vertex_rotate_seq.resize(m+n);

    for (int i=0; i<n; ++i) {
      vertex_coord_seq[i+m].resize(vertex_num);
      vertex_rotate_seq[i+m].resize(vertex_num);
      for (int j=0; j<vertex_num; ++j) {
	float x,y,z;
	fid >> x >> y >> z;
	vertex_coord_seq[i+m][j]=Point(x,y,z);	
      }
    }

    fid.close();
      
  }

  void ManifoldTriMesh::load_examples(std::string filename) {
    filename.append(".rot");
    std::fstream fid; fid.open(filename.c_str());
    int vn, m; fid >> vn >> m;
    vertex_rotate_exp.resize(m);
    for (int i = 0; i < m; ++i) vertex_rotate_exp[i].resize(vn);

    for (int i = 0; i < vn; ++i) {
      for (int j = 0; j < m; ++j) {
	float x,y,z;
	fid >> x >> y >> z;
	vertex_rotate_exp[j][i] = Vector(x,y,z);
      }	
    }

    fid.close();
  }

  inline void assembly_rotation_matrix(float * const &A, const float &w, const Point &a, const Point &b, const Vector &v) {
    A[0] += w * (a.x() - b.x()) * v.x();
    A[1] += w * (a.x() - b.x()) * v.y();
    A[2] += w * (a.x() - b.x()) * v.z();
    A[3] += w * (a.y() - b.y()) * v.x();
    A[4] += w * (a.y() - b.y()) * v.y();
    A[5] += w * (a.y() - b.y()) * v.z();
    A[6] += w * (a.z() - b.z()) * v.x();
    A[7] += w * (a.z() - b.z()) * v.y();
    A[8] += w * (a.z() - b.z()) * v.z();

  }


  inline Vector solve_rotation(float * const &A) {//robust implement computing axis-angle rotation
    float U[9], V[9], S[3];
    float norm=0;
    for (int i=0;i<9;++i) norm+=A[i]*A[i]; norm = 1/std::sqrt(norm);
    for (int i=0;i<9;++i) A[i]*=norm;
    
    fastsvd(A, U, V, S);

    A[0] = V[0]*U[0] + V[3]*U[3] + V[6]*U[6];
    A[1] = V[1]*U[0] + V[4]*U[3] + V[7]*U[6];
    A[2] = V[2]*U[0] + V[5]*U[3] + V[8]*U[6];
    A[3] = V[0]*U[1] + V[3]*U[4] + V[6]*U[7];
    A[4] = V[1]*U[1] + V[4]*U[4] + V[7]*U[7];
    A[5] = V[2]*U[1] + V[5]*U[4] + V[8]*U[7];
    A[6] = V[0]*U[2] + V[3]*U[5] + V[6]*U[8];
    A[7] = V[1]*U[2] + V[4]*U[5] + V[7]*U[8];
    A[8] = V[2]*U[2] + V[5]*U[5] + V[8]*U[8];

    float trace = A[0] + A[4] + A[8], theta = std::acos((trace-1)/2.);
    const float epsilon = 1E-6;
    
    if (trace >= 3 -epsilon) 
      return  (.5 - (trace -3)/12.) * Vector(A[5]-A[7], A[6]-A[2], A[1]-A[3]);
    else if (trace > -1 + epsilon)
      return  (.5*theta/std::sin(theta)) * Vector(A[5]-A[7], A[6]-A[2], A[1]-A[3]);
    else  if (A[0]>A[4] && A[0] > A[8]) {
      float s = std::sqrt(A[0]-A[4]-A[8]+1);
      Vector v = Vector(s, (A[1]+A[3])/s, (A[2]+A[6])/s);
      norm = std::sqrt(v*v);
      return MESHTK_PI * v/norm;
    } else if (A[4]> A[8]) {
      float s = std::sqrt(-A[0]+A[4]-A[8]+1);
      Vector v = Vector((A[1]+A[3])/s, s, (A[5]+A[7])/s);
      norm = std::sqrt(v*v);
      return MESHTK_PI * v/norm;
    } else {
      float s = std::sqrt(-A[0]-A[4]+A[8]+1);
      Vector v = Vector((A[2]+A[6])/s, (A[5]+A[7])/s, s);
      norm = std::sqrt(v*v);
      return MESHTK_PI * v/norm;      
    }

  }

  void ManifoldTriMesh::compute_rotate_sequence() {
    int n = vertex_rotate_seq.size();
    
    for (int j=0; j < vertex_num; ++j) {
      for (int i=0; i < n; ++i) {
	float A[9]={};
	HV_circulator hv = IV[j]->vertex_begin();
	do {
	  int id, ida, idb; float w; Halfedge_handle h;

	  h = hv;
	  id = h->index; 
	  ida = h->vertex()->index;
	  idb = h->prev()->vertex()->index;
	  w = std::fabs(halfedge_cot[id]); 
	  assembly_rotation_matrix(A, w, vertex_coord_seq[i][ida], vertex_coord_seq[i][idb], halfedge_vec[id]);

	  h = hv->prev();
	  id = h->index; 
	  ida = h->vertex()->index;
	  idb = h->prev()->vertex()->index;
	  w = std::fabs(halfedge_cot[id]); 
	  assembly_rotation_matrix(A, w, vertex_coord_seq[i][ida], vertex_coord_seq[i][idb], halfedge_vec[id]);	  

	  h = hv->next();
	  id = h->index; 
	  ida = h->vertex()->index;
	  idb = h->prev()->vertex()->index;
	  w = std::fabs(halfedge_cot[id]); 
	  assembly_rotation_matrix(A, w, vertex_coord_seq[i][ida], vertex_coord_seq[i][idb], halfedge_vec[id]);	  
	  
	  ++hv;
	} while (hv != IV[j]->vertex_begin());
	vertex_rotate_seq[i][j] = solve_rotation(A);
      }
    }

    for (int i=0; i < n; ++i) {
      BooleanFunction check; check.resize(vertex_num);
      std::vector<int> label; label.resize(vertex_num);

      int head, tail=0;
      do {
	int iter = 0; while (check[iter]) ++iter;
	check[iter] = true;//reference
	
	label[head = tail++] = iter;
	while (head != tail) {
	  Vector vec_head = vertex_rotate_seq[i][head], vec_head_normalized;
	  double norm_head = std::sqrt(vec_head*vec_head);
	  if (norm_head > 1E-6) vec_head_normalized = vec_head/norm_head;

	  HV_circulator hv = IV[label[head]]->vertex_begin();
	  do {
	    int idx = hv->prev()->vertex()->index;
	    if (!check[idx]) {
	      Vector v = vertex_rotate_seq[i][idx];
	      double norm_tail = std::sqrt(v*v);
	      int half_step=0;
	      if (norm_tail < 1E-6) {
		if (norm_head > 1E-6) {		  
		  double vec_tail_project = v * vec_head_normalized;
		  half_step = int ((norm_head - vec_tail_project)/MESHTK_PI);
		  half_step = std::floor((half_step+1.)/2.);		  
		  if (half_step) vertex_rotate_seq[i][idx] = vertex_rotate_seq[i][idx] + half_step * 2* MESHTK_PI * vec_head_normalized;
		}
	      } else {
		v = v/norm_tail;
		double vec_head_project =  vec_head * v;
		half_step = int ((vec_head_project - norm_tail)/MESHTK_PI);
		half_step = std::floor((half_step +1.)/2.);
		if (half_step) vertex_rotate_seq[i][idx] = vertex_rotate_seq[i][idx] + half_step * 2* MESHTK_PI * v;
	      }
	      
	      if (half_step) std::cout << "+";
	      check[idx] = true;
	      label[tail++] = idx;
	    }
	    ++hv;
	  } while (hv != IV[label[head]]->vertex_begin());
	  ++head;
	}
      } while (head < vertex_num);

    }

  }

  void ManifoldTriMesh::print_rotate_sequence(std::string filename) {
    int m = vertex_rotate_seq.size();
    filename.append(".rot");

    std::ofstream fid; fid.open(filename.c_str()); fid << vertex_num << "\t" << 3*m << std::endl;
    for (int j=0; j<vertex_num; ++j) {
      for (int i=0; i<m; ++i) 
	fid << vertex_rotate_seq[i][j].x() << "\t" 
	    << vertex_rotate_seq[i][j].y() << "\t"
	    << vertex_rotate_seq[i][j].z() << "\t";
      fid << std::endl;
    }

  }

  void ManifoldTriMesh::load_proxy_bone(std::string filename) {
    std::string colorfile=filename;
    colorfile.append(".color");
    filename.append(".id");
    std::fstream fid; fid.open(filename.c_str());
    vertex_label.resize(vertex_num);
    for (int i=0; i<vertex_num; ++i) fid >> vertex_label[i];
    fid.close();    

    fid.open(colorfile.c_str());
    float c;
    while (fid >> c) label_color.push_back(c);
    fid.close();
  }

}
