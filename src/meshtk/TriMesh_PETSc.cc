/*
  FILE: meshtk/TriMesh_PETSc.cc This file is part of tri-mesh-toolkit.
  It is a C++ header file which implement PETSc matrix for mesh.
  
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

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"
#include <petscmat.h>


#include <iostream>
#include <fstream>

namespace meshtk {

  // mass matrix and stiff matrix of Laplace Beltrami operator
  static Mat mass_mat, stiff_mat;

  static PetscInt mat_size;

  // eigenvalues and eigenvectors of generalized eigen problem defined by FEM_LBmat
  static std::vector<double> eig_value;
  static std::vector<Vec> eig_vector;
  static std::vector<double> eig_vector_sqr_norm;
  // number of eigenvalues computed
  static PetscInt eig_num;



  static char help[] = "";

  void TriMesh::PETSc_init(int argc, char **argv){
    PetscInitialize(&argc,&argv,(char *)0,help);
  }

  void TriMesh::PETSc_destroy() {
    
    MatDestroy(&mass_mat); MatDestroy(&stiff_mat);    
    
    for (int i = 0; i < eig_num; ++i) VecDestroy(&eig_vector[i]);
    
  }



const PetscInt I1[100]={
	0,	7,	-7,	57,	-24,	0,	0,	24,	-57,	0,
	7,	0,	-7,	-24,	57,	-57,	24,	0,	0,	0,
	-7,	-7,	68,	-6,	-6,	30,	-51,	-51,	30,	0,
	57,	-24,	-6,	135,	54,	27,	27,	27,	-135,	-162,
	-24,	57,	-6,	54,	135,	-135,	27,	27,	27,	-162,
	0,	-57,	30,	27,	-135,	135,	-108,	-27,	-27,	162,
	0,	24,	-51,	27,	27,	-108,	135,	135,	-27,	-162,
	24,	0,	-51,	27,	27,	-27,	135,	135,	-108,	-162,
	-57,	0,	30,	-135,	27,	-27,	-27,	-108,	135,	162,
	0,	0,	0,	-162,	-162,	162,	-162,	-162,	162,	324};


const PetscInt I2[100]={
	0,	-7,	7,	-57,	24,	0,	0,	-24,	57,	0,
	-7,	68,	-7,	30,	-51,	-51,	30,	-6,	-6,	0,
	7,	-7,	0,	0,	0,	24,	-57,	57,	-24,	0,
	-57,	30,	0,	135,	-108,	-27,	-27,	27,	-135,	162,
	24,	-51,	0,	-108,	135,	135,	-27,	27,	27,	-162,
	0,	-51,	24,	-27,	135,	135,	-108,	27,	27,	-162,
	0,	30,	-57,	-27,	-27,	-108,	135,	-135,	27,	162,
	-24,	-6,	57,	27,	27,	27,	-135,	135,	54,	-162,
	57,	-6,	-24,	-135,	27,	27,	27,	54,	135,	-162,
	0,	0,	0,	162,	-162,	-162,	162,	-162,	-162,	324};

const PetscInt I3[100]={
	68,	-7,	-7,	-51,	30,	-6,	-6,	30,	-51,	0,
	-7,	0,	7,	24,	-57,	57,	-24,	0,	0,	0,
	-7,	7,	0,	0,	0,	-24,	57,	-57,	24,	0,
	-51,	24,	0,	135,	-108,	27,	27,	-27,	135,	-162,
	30,	-57,	0,	-108,	135,	-135,	27,	-27,	-27,	162,
	-6,	57,	-24,	27,	-135,	135,	54,	27,	27,	-162,
	-6,	-24,	57,	27,	27,	54,	135,	-135,	27,	-162,
	30,	0,	-57,	-27,	-27,	27,	-135,	135,	-108,	162,
	-51,	0,	24,	135,	-27,	27,	27,	-108,	135,	-162,
	0,	0,	0,	-162,	162,	-162,	-162,	162,	-162,	324};

const PetscInt I4[100]
  ={  76,  11,  11,  18,   0,  27,  27,   0,  18,  36,
      11,  76,  11,   0,  18,  18,   0,  27,  27,  36,
      11,  11,  76,  27,  27,   0,  18,  18,   0,  36,
      18,   0,  27, 540,-189,-135, -54,-135, 270, 162,
       0,  18,  27,-189, 540, 270,-135, -54,-135, 162,
      27,  18,   0,-135, 270, 540,-189,-135, -54, 162,
      27,   0,  18, -54,-135,-189, 540, 270,-135, 162,
       0,  27,  18,-135, -54,-135, 270, 540,-189, 162,
      18,  27,   0, 270,-135, -54,-135,-189, 540, 162,
      36,  36,  36, 162, 162, 162, 162, 162, 162,1944};  




  void TriMesh::PETSc_assemble_cubicFEM_LBmat(){
    clock_start("Assemble cubic FE");

    mat_size = vertex_num + halfedge_num + facet_num;

    PetscInt *nnz =new PetscInt[mat_size];

    for (PetscInt i=0;i<vertex_num;++i)
      nnz[i]=IV[i]->vertex_degree()*6+1;

    for (PetscInt i=0;i<halfedge_num;i++)
      nnz[i+vertex_num]=12;

    for (PetscInt i=0;i<facet_num;i++)
      nnz[i+vertex_num+halfedge_num]=1;


    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &mass_mat);
    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &stiff_mat);

    delete [] nnz;

    MatSetOption(mass_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
    MatSetOption(stiff_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);


    for (PetscInt i = 0; i < facet_num; ++i) {
      Halfedge_handle h = IF[i]->halfedge();
      
      double l1 = halfedge_length[h->next()->index],
	l2 = halfedge_length[h->index],
	l3 = halfedge_length[h->prev()->index];

      int idx[10];
      idx[0]=h->vertex()->index;
      idx[1]=h->next()->vertex()->index;
      idx[2]=h->prev()->vertex()->index;
      idx[3]=h->next()->index + vertex_num;
      idx[4]=h->next()->opposite()->index + vertex_num;
      idx[5]=h->prev()->index + vertex_num;
      idx[6]=h->prev()->opposite()->index + vertex_num;
      idx[7]=h->index + vertex_num;
      idx[8]=h->opposite()->index + vertex_num;
      idx[9]=i + vertex_num + halfedge_num;

      double mass[100], stiff[100];
      for (PetscInt j = 0; j < 100; ++j) {
	mass[j]=facet_area[i]*I4[j]/6720.;
	stiff[j]=(l1*l1*I1[j]+l2*l2*I2[j]+l3*l3*I3[j])/(320.* facet_area[i]);
      }

      MatSetValues(stiff_mat, 10, idx, 10, idx, stiff, ADD_VALUES);
      MatSetValues(mass_mat, 10, idx, 10, idx, mass, ADD_VALUES);	       
      
    }


    MatAssemblyBegin(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mass_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mass_mat, MAT_FINAL_ASSEMBLY);



    clock_end();
  }

  void TriMesh::PETSc_assemble_linearFEM_LBmat(){

  }
  
  void PETSc_export_mat(const char *filename, Mat mat) {
    PetscViewer fviewer;
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &fviewer);
    MatView(mat, fviewer);
    PetscViewerDestroy(&fviewer);
  }

  void TriMesh::PETSc_load_LBmat(std::string name) {
    
    PetscViewer fmass, fstiff;

    std::string filename_stiff = name, filename_mass = name;
    filename_stiff += ".stiff"; filename_mass += ".mass";

    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename_mass.c_str() , FILE_MODE_READ, &fmass);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename_stiff.c_str(), FILE_MODE_READ, &fstiff);

    MatCreate(PETSC_COMM_SELF, &mass_mat);
    MatCreate(PETSC_COMM_SELF, &stiff_mat);
    MatSetType(mass_mat, MATSEQSBAIJ);
    MatSetType(stiff_mat, MATSEQSBAIJ);

    MatLoad(mass_mat,  fmass);
    MatLoad(stiff_mat, fstiff);
      
    PetscViewerDestroy(&fmass); PetscViewerDestroy(&fstiff);
    
  }

  void TriMesh::PETSc_export_LBmat( std::string name) {

    std::string filename_stiff = name, filename_mass = name;
    filename_stiff += ".stiff"; filename_mass += ".mass";
    PETSc_export_mat(filename_mass.c_str(), mass_mat);
    PETSc_export_mat(filename_stiff.c_str(), stiff_mat);

  }

  

  PetscScalar vec_inner_prod(Vec &v1, Vec &v2) {
    Vec tmp;
    PetscScalar inner_prod;
    VecDuplicate(v1, &tmp);

    MatMult(mass_mat, v2, tmp);
    VecDot(v1, tmp, &inner_prod);
    return inner_prod;
  }

  void TriMesh::PETSc_load_LBeigen(std::string name) {

    clock_start("Load eigen");

    std::ifstream f_eigvalue;
    double eig;
    std::string eigvalue_filename = name;
    eigvalue_filename.append("/_ev.ascii");

    f_eigvalue.open(eigvalue_filename.c_str());
    
    if (!f_eigvalue.is_open()) {
      std::cerr << "Error: Cannot open file " << eigvalue_filename << std::endl;
      exit(1);
    }
    
    eig_value.clear();
    while (!(f_eigvalue >> eig).eof()) { eig_value.push_back(eig); } 
    f_eigvalue.close();
    eig_num = eig_value.size();
    
    eig_vector.resize(eig_num);
    eig_vector_sqr_norm.resize(eig_num);

    for (int i =0; i< eig_num; ++i) {
      char ffile[PETSC_MAX_PATH_LEN];
      PetscViewer fd;
  
      sprintf(ffile,"%s/_%d.petsc",name.c_str(),i);
      PetscViewerBinaryOpen(PETSC_COMM_SELF, ffile, FILE_MODE_READ, &fd);
      VecCreate(PETSC_COMM_SELF, &eig_vector[i]);
      VecSetType(eig_vector[i], VECSEQ);
      VecLoad(eig_vector[i], fd);
      PetscViewerDestroy(&fd);
      
      eig_vector_sqr_norm[i] = vec_inner_prod(eig_vector[i], eig_vector[i]);
    }

    clock_end();

  }

  double TriMesh::update_vertex_biharmonic(int source_vertex_index,
					   ScalarFunction & biharmonic_distance) {
 

    for (int j=0; j < vertex_num; ++j) biharmonic_distance[j] =0;

    for (int i=eig_num-1; i> 0; --i) {
      PetscScalar *vector;
      VecGetArray(eig_vector[i], &vector);
      
      PetscScalar target = vector[source_vertex_index], sqr_eigenvalue = eig_value[i] * eig_value[i]; 
      for (int j=0; j < vertex_num; ++j) {
	double tmp = vector[j] -target;
	biharmonic_distance[j] += tmp * tmp / eig_vector_sqr_norm[i] / sqr_eigenvalue;
      }
    }

    double total_area_distance=0;
    for (int j=0; j < vertex_num; ++j) {

      biharmonic_distance[j] = std::sqrt(biharmonic_distance[j]);
      total_area_distance += biharmonic_distance[j] * vertex_area[j];
    }

    
    return total_area_distance / total_area;
  }

  

  int TriMesh::PETSc_assemble_export_BiHDM(std::vector<int> & sampling, //init sampling provided, keypoint based
				     int addition_size, // expect additional sampling size
				     double stop_criterion) {    
    clock_start("Assemble Nystrom sampling of BiHDM");


    int init_size = sampling.size();
    int total_size = init_size;
    
    //double max_A_distance;
    double max_B_distance = 0;
    int max_B_distance_index =  vertex_num * rand()/ RAND_MAX;

    ScalarFunction distance;
    ScalarFunction new_distance;
    
    distance.resize(vertex_num);



    for (int i =0;i<vertex_num; ++i) distance[i] =-1;


    for (int i = 0; i < init_size; ++i ){
      new_distance.clear();
      new_distance.resize(vertex_num);
      update_vertex_biharmonic(sampling[i], new_distance);
      
      for (int j=0; j< vertex_num; ++j) 
	if (new_distance[j] < distance[j] || distance[j] <0) {
	  distance[j] = new_distance[j];
	}
      //export  sampling[i] and new_distance here
    }

    for (int j=0; j<vertex_num; ++j) 
      if (distance[j] > max_B_distance) {
	max_B_distance = distance[j];
	max_B_distance_index = j;
      }

    

    while (total_size < (init_size + addition_size)){
      new_distance.clear();
      new_distance.resize(vertex_num); 
      sampling.push_back(max_B_distance_index);
      update_vertex_biharmonic(max_B_distance_index, new_distance);

      for (int j=0; j< vertex_num; ++j) 
	if (new_distance[j] < distance[j] || distance[j] <0) {
	  distance[j] = new_distance[j];
	}
      //export max_B_distance_index and new_distance here

      max_B_distance = 0;
      for (int j=0; j<vertex_num; ++j) 
	if (distance[j] > max_B_distance) {
	  max_B_distance = distance[j];
	  max_B_distance_index = j;
	}
      ++ total_size;

    }

    clock_end();


    return total_size;
  }

  void TriMesh::PETSc_assemble_Fourier_BiHDM() {
  }
  void TriMesh::PETSc_export_Fourier_BiHDM() {
  }


}




