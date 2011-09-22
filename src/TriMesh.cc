
#include <iostream>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"

using namespace meshtk;

void TriMesh::read(std::string file, std::string type){
  if (type.compare("off")==0){
    std::ifstream mesh_Fin;
    file.append("."); file.append("off");    
    mesh_Fin.open(file.c_str());
    if (!mesh_Fin.is_open())
      std::cerr << "Cannot open file " << file << std::endl;
    mesh_Fin >> P;// mesh read
    mesh_Fin.close();
  }


  if (!P.is_pure_triangle()) {std::cerr<< "Error: input mesh has non-triangle facet!" <<std::endl; exit(1);}
  if (!P.is_closed()) {std::cout<< "Warning: input mesh seems not to be watertight" <<std::endl;} 
};

void TriMesh::write(std::string file, std::string type){
  if (type.compare("off")==0){
    std::ofstream mesh_Fout;
    file.append("."); file.append(type);
    mesh_Fout.open(file.c_str());
    mesh_Fout << P;// mesh output
    std::cout << "Export mesh to: " << file <<  std::endl;
    mesh_Fout.close();
  }
};

void TriMesh::init_index(){
  //int n=(P.size_of_halfedges()+P.size_of_border_edges())/2;

  halfedge_num = P.size_of_halfedges();
  vertex_num = P.size_of_vertices();
  facet_num = P.size_of_facets();

  IH = ISHalfedge_list(halfedge_num);
  IV =  ISVertex_list(vertex_num);
  IF =  ISFacet_list(facet_num);

  halfedge_vec.resize(halfedge_num);
  vertex_norm.resize(vertex_num);
  facet_norm.resize(facet_num);
  vertex_area.resize(vertex_num);
  facet_area.resize(facet_num);

  vertex_avg_len.resize(vertex_num);

  int index_count=0;
  for(Vertex_iterator vitr= P.vertices_begin();vitr!= P.vertices_end();
      IV[index_count]=vitr, vitr->index = index_count++, vitr++);
  index_count=0;
  for(Halfedge_iterator eitr= P.halfedges_begin();eitr!= P.halfedges_end();
      IH[index_count]=eitr, eitr->index = index_count++, eitr++);
  index_count=0;
  for(Facet_iterator fitr= P.facets_begin(); fitr!= P.facets_end(); 
      IF[index_count]=fitr, fitr->index = index_count++, fitr++);
};


double TriMesh::update_halfedge(){

  double avg_len=0;
  Halfedge_handle h;
  for (int i=0;i<halfedge_num;i++) 
    {
      h = IH[i];
      halfedge_vec[i] = h->vertex()->point() - h->prev()->vertex()->point();
      avg_len += CGAL::sqrt(halfedge_vec[i] * halfedge_vec[i]);
    }
  return avg_edge_len = avg_len/halfedge_num;
}


double TriMesh::update_facet(){
  Vector normal;
  Halfedge_handle h;

  total_area=0;
  for (int i=0;i<facet_num;i++){
    h=IF[i]->halfedge();
    normal = CGAL::cross_product(halfedge_vec[h->index], halfedge_vec[h->next()->index]);
    total_area += (facet_area[i] = normal * normal) /2. ;
    facet_norm[i] = normal / CGAL::sqrt(facet_area[i]); 
    facet_area[i] /= 2.;
  }
  
  return  total_area;
}


void TriMesh::update_vertex(){
  /*
  double sigma = 2 * avg_edge_len;

  for (int i=0;i<vertex_num;i++){
    Vector normal(0,0,0), tmp;
    double scale;
    HV_circulator hv=IV[i]->vertex_begin();
    do{
      tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
      scale = CGAL::sqrt(tmp * tmp);
      normal = normal + std::exp(- scale * scale / (sigma *sigma)) * facet_norm[hv->facet()->index];
    }while (++hv!=IV[i]->vertex_begin());
    vertex_norm[i] = normal / CGAL::sqrt(normal * normal);
  }
  */


  for (int i=0;i<vertex_num;i++) {
    Point p= IV[i]->point();
    if (p.x() < coord_min_x || i==0) coord_min_x = p.x();
    if (p.y() < coord_min_y || i==0) coord_min_y = p.y();
    if (p.z() < coord_min_z || i==0) coord_min_z = p.z();
    if (p.x() > coord_max_x || i==0) coord_max_x = p.x();
    if (p.y() > coord_max_y || i==0) coord_max_y = p.y();
    if (p.z() > coord_max_z || i==0) coord_max_z = p.z();
    
  }

  for (int i=0; i<vertex_num; i++){
    HV_circulator hv=IV[i]->vertex_begin();
    double area = 0, total_len =0;
    Vector tmp;
    int k=0;

    do {
      if (hv->facet()==NULL) continue;
      area += facet_area[hv->facet()->index];

      tmp = halfedge_vec[hv->index];
      total_len += CGAL::sqrt(tmp * tmp);		    
      ++k;      

    }while (++hv!=IV[i]->vertex_begin());
    vertex_area[i] = area/3.;
    vertex_avg_len[i] = total_len / k;

  }

  facet2vertex_point_average( facet_norm, vertex_norm, Vector(0,0,0));
  for (int i=0;i<vertex_num;i++){
    vertex_norm[i] = vertex_norm[i] / CGAL::sqrt(vertex_norm[i] * vertex_norm[i]);
  }

    
}


void TriMesh::update_base(){//base update halfedge, facet, vertex.

  update_halfedge();

  update_facet();
  
  update_vertex();

};


void TriMesh::update_vertex_localchart(){
  vertex_LC[0] = Vec_Fun(vertex_num);
  vertex_LC[1] = Vec_Fun(vertex_num);
  
  for (int i=0;i<vertex_num;i++) localchart(vertex_LC[0][i],vertex_LC[1][i], vertex_norm[i]);
}

void TriMesh::update_facet_localchart(){
  facet_LC[0] = Vec_Fun(facet_num);
  facet_LC[1] = Vec_Fun(facet_num);

  for (int i=0;i<facet_num;i++) localchart(facet_LC[0][i],facet_LC[1][i], facet_norm[i]);
}


Curvature prin_curv(double e, double f, double g){
  double H= (e+g)/2;
  double G=std::sqrt(4*f*f+(e-g)*(e-g))/2;

  return meshtk::Curvature(H+G, H-G, H, H*H-G*G);
}

void TriMesh::update_facet_curvature(){

  facet_CT[0].resize(facet_num);
  facet_CT[1].resize(facet_num);
  facet_CT[2].resize(facet_num);

  facet_PC0.resize(facet_num);
  facet_PC1.resize(facet_num);
  facet_hcurv.resize(facet_num);
  facet_kcurv.resize(facet_num);
  
  //for (int i=0;i<1;i++){
  for (int i=0;i<facet_num;i++){

    //compute curvature for each face

    //step 1. assembly matrix
    Halfedge_handle h[3];
    h[0]=IF[i]->halfedge();
    h[1]=h[0]->next();
    h[2]=h[1]->next();

    Vector e[3],n[3];
    for (int j=0;j<3;j++) e[j] = halfedge_vec[h[j]->index];
    
    for (int j=0;j<3;j++) n[j] = vertex_norm[h[j]->vertex()->index] - vertex_norm[h[j]->prev()->vertex()->index];

    double e_LC[3][2], n_LC[3][2];
    for (int j=0;j<3;j++) {
      localcoord(e[j],facet_LC[0][i],facet_LC[1][i], e_LC[j]);
      localcoord(n[j],facet_LC[0][i],facet_LC[1][i], n_LC[j]);
    }


    double matrix_l[2][2]={}, matrix_r[2][2]={};
    for (int j=0;j<2;j++)
      for (int k=0;k<2;k++)
	for (int l=0;l<3;l++){
	  matrix_l[j][k]+=e_LC[l][j]*e_LC[l][k];
	  matrix_r[j][k]+=n_LC[l][j]*e_LC[l][k];
	}
    
    double final_linv[3][3]={}, final_r[3];
    double a,b,c, deter;
    a = matrix_l[0][0]; b = matrix_l[0][1]; c = matrix_l[1][1];
    deter=(a+c)*(a*c-b*b);
    
    final_r[0] = matrix_r[0][0];
    final_r[1] = matrix_r[0][1]+matrix_r[1][0];
    final_r[2] = matrix_r[1][1];

    
    final_linv[0][0]=(a+c)*c-b*b;
    final_linv[0][1]=final_linv[1][0]=-b*c;
    final_linv[0][2]=final_linv[2][0]=b*b;
    final_linv[1][1]=a*c;
    final_linv[1][2]=final_linv[2][1]=-b*a;
    final_linv[2][2]=(a+c)*a-b*b;

    
    facet_CT[0][i] = (final_linv[0][0]*final_r[0] + final_linv[0][1]*final_r[1] + final_linv[0][2]*final_r[2])/deter;
    
    facet_CT[1][i] = (final_linv[1][0]*final_r[0] + final_linv[1][1]*final_r[1] + final_linv[1][2]*final_r[2])/deter;
    facet_CT[2][i] = (final_linv[2][0]*final_r[0] + final_linv[2][1]*final_r[1] + final_linv[2][2]*final_r[2])/deter;        


    Curvature curv = prin_curv(facet_CT[0][i], facet_CT[1][i], facet_CT[2][i]);

    facet_PC0[i] = curv.PC0;
    facet_PC1[i] = curv.PC1;
    facet_hcurv[i] = curv.hcurv;
    facet_kcurv[i] = curv.kcurv;


    //step 2. solve matrix

  }
}



void TriMesh::update_vertex_curvature(){
  vertex_CT[0].resize(vertex_num);
  vertex_CT[1].resize(vertex_num);
  vertex_CT[2].resize(vertex_num);

  vertex_PC0.resize(vertex_num);
  vertex_PC1.resize(vertex_num);

  vertex_hcurv.resize(vertex_num);
  vertex_kcurv.resize(vertex_num);
  double sigma = avg_edge_len;
  Vector tmp;

  for (int i=0; i < vertex_num; ++i) {
    HV_circulator hv = IV[i]->vertex_begin();
    double total_scale=0, scale;
    vertex_CT[0][i] = vertex_CT[1][i] = vertex_CT[2][i] =0;

    sigma = vertex_avg_len[i];

    do{
      if (hv->facet() == NULL) continue;

      scale = std::exp(- scale * scale / (2 * sigma * sigma));
      total_scale += scale;

      double x[2][2]; 
      int j=hv->facet()->index;

      for (int m=0;m<2;++m) for (int k=0;k<2;++k) x[m][k] = vertex_LC[m][i] * facet_LC[k][j];

      
      vertex_CT[0][i] += scale * 
	(facet_CT[0][j]*x[0][0]*x[0][0] + 2*facet_CT[1][j]*x[0][0]*x[0][1] + facet_CT[2][j]*x[0][1]*x[0][1]);
      vertex_CT[1][i] += scale *
	(facet_CT[0][j]*x[0][0]*x[1][0] + facet_CT[1][j]*(x[0][0]*x[1][1] + x[0][1]*x[1][0])+ facet_CT[2][j]*x[0][1]*x[1][1]);
      vertex_CT[2][i] += scale *
	(facet_CT[0][j]*x[1][0]*x[1][0] + 2*facet_CT[1][j]*x[1][0]*x[1][1] + facet_CT[2][j]*x[1][1]*x[1][1]);

    }while (++hv != IV[i]->vertex_begin());

    vertex_CT[0][i]/=total_scale;
    vertex_CT[1][i]/=total_scale;
    vertex_CT[2][i]/=total_scale;

 
    Curvature curv = prin_curv(vertex_CT[0][i], vertex_CT[1][i], vertex_CT[2][i]);

    vertex_PC0[i] = curv.PC0;
    vertex_PC1[i] = curv.PC1;
    vertex_hcurv[i] = curv.hcurv;
    vertex_kcurv[i] = curv.kcurv;
  }

  /*
  facet2vertex_point_average( facet_hcurv, vertex_hcurv, 0.);
  facet2vertex_point_average( facet_kcurv, vertex_kcurv, 0.);
  */
}

void TriMesh::update_curvature(){

  update_facet_localchart();
  update_facet_curvature();

  update_vertex_localchart();
  update_vertex_curvature();

}

