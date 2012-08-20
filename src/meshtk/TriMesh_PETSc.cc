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

#include <assert.h>

namespace meshtk {

  static PetscInt FEM_degree = 3;//cubic FEM in default;

  // index map
  static PetscInt* FEM_index_map;

  // mass matrix and stiff matrix of Laplace Beltrami operator
  static Mat mass_mat, stiff_mat;
  static PetscInt mat_size;


  // Fourier bent bihamonic matrix 
  static Mat fbihd_mat;
  static Vec bihd_trace;

  // eigenvalues and eigenvectors of generalized eigen problem defined by FEM_LBmat
  static std::vector<double> eig_value;
  static std::vector<Vec> eig_vector;
  static std::vector<double> eig_vector_sqr_norm;
  // number of eigenvalues computed
  static PetscInt eig_num;

  //static std::vector<std::vector<double> > vertex_hks;
  //static int vertex_hks_dim;







  static char help[] = "";

  void TriMesh::PETSc_init(int argc, char **argv){
    PetscInitialize(&argc,&argv,(char *)0,help);
    

  }

  void TriMesh::PETSc_destroy() {
    
    MatDestroy(&mass_mat); MatDestroy(&stiff_mat);    
    
    for (int i = 0; i < eig_num; ++i) VecDestroy(&eig_vector[i]);

    MatDestroy(&fbihd_mat);
    VecDestroy(&bihd_trace);

    delete [] FEM_index_map;
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


  int TriMesh::map_cubicFEM_indices(bool dirichlet) {
    //only to update FEM_index_map array    
    FEM_index_map = new PetscInt[vertex_num + halfedge_num + facet_num];
    for (PetscInt i=0; i<vertex_num + halfedge_num + facet_num; ++i) FEM_index_map[i] = 0;

    if (dirichlet) {
      for (int i=0;i<halfedge_num;i++) 
	{
	  Halfedge_handle h = IH[i];
	  if (h->is_border())
	    {
	      FEM_index_map[i+vertex_num] = -1;
	      FEM_index_map[h->opposite()->index+vertex_num] = -1;
	      FEM_index_map[h->vertex()->index] = -1;
	      FEM_index_map[h->opposite()->facet()->index + vertex_num + halfedge_num] =-1;
	    }
	}
    }
    
    int count =0;
    int vertex_count = 0, halfedge_count=0, facet_count =0;
    for (PetscInt i=0; i<vertex_num; ++i)
      if (FEM_index_map[i]==0)
	FEM_index_map[i] = count++;
    vertex_count = count;

    for (PetscInt i=vertex_num; i<vertex_num + halfedge_num; ++i)
      if (FEM_index_map[i]==0)
	FEM_index_map[i] = count++;
    halfedge_count = count - vertex_count;

    for (PetscInt i=vertex_num+halfedge_num; i<vertex_num + halfedge_num + facet_num; ++i)
      if (FEM_index_map[i]==0)
	FEM_index_map[i] = count++;
    facet_count = count - vertex_count - halfedge_count;

     return count; 

  }

  void TriMesh::PETSc_assemble_cubicFEM_LBmat(bool dirichlet){
    clock_start("Assemble cubic FE");


    // FEM index map
    mat_size = map_cubicFEM_indices(dirichlet);
    
    PetscInt *nnz =new PetscInt[mat_size];

    for (PetscInt i=0;i<vertex_num;++i)
      if (FEM_index_map[i]>=0) nnz[FEM_index_map[i]]=IV[i]->vertex_degree()*6+1;

    for (PetscInt i=vertex_num;i<vertex_num + halfedge_num;i++)
      if (FEM_index_map[i]>=0) nnz[FEM_index_map[i]]=12;

    for (PetscInt i=vertex_num + halfedge_num;i<vertex_num + halfedge_num + facet_num;i++)
      if (FEM_index_map[i]>=0) nnz[FEM_index_map[i]]=1;


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

      for (PetscInt j=0; j < 10; j++) idx[j] = FEM_index_map[idx[j]];

      PetscScalar mass[100], stiff[100];
      for (PetscInt j = 0; j < 100; ++j) {
	mass[j]= facet_area[i]*I4[j]/6720.;
	stiff[j]= (l1*l1*I1[j]+l2*l2*I2[j]+l3*l3*I3[j])/(320.* facet_area[i]);
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


const PetscInt J1[9]
= {2, -1, -1, 
   -1, 0, 1, 
   -1, 1, 0};
const PetscInt J2[9]
= {0, -1, 1, 
   -1, 2, -1, 
   1, -1, 0};
const PetscInt J3[9]
= {0, 1, -1,
   1, 0, -1, 
   -1, -1, 2};
const PetscInt J4[9]
= {2, 1, 1, 1, 2, 1, 1, 1, 2};

  int TriMesh::map_linearFEM_indices(bool dirichlet) {
    //only to update FEM_index_map array
    FEM_index_map = new PetscInt[vertex_num];
    for (PetscInt i=0; i<vertex_num; ++i) FEM_index_map[i] = 0;

    if (dirichlet) {
      for (int i=0;i<halfedge_num;i++) 
	{
	  Halfedge_handle h = IH[i];
	  if (h->is_border())
	    {
	      FEM_index_map[h->vertex()->index] = -1;
	    }
	}
    }
    
    int count =0;
    for (PetscInt i=0; i<vertex_num; ++i)
      if (FEM_index_map[i]==0)
	FEM_index_map[i] = count++;

     return count; 
  }


  void TriMesh::PETSc_assemble_linearFEM_LBmat(bool dirichlet){
    clock_start("Assemble linear FE");

    mat_size = map_linearFEM_indices(dirichlet);
    

    PetscInt *nnz =new PetscInt[mat_size];
    for (PetscInt i=0;i<vertex_num;++i)
      if (FEM_index_map[i] >=0) nnz[FEM_index_map[i]]=IV[i]->vertex_degree()+1;  

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &mass_mat);
    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &stiff_mat);

    delete [] nnz;
    MatSetOption(mass_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
    MatSetOption(stiff_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    for (PetscInt i = 0; i < facet_num; ++i) {
      Halfedge_handle h = IF[i]->halfedge();
      
      double l1 = halfedge_length[h->prev()->index],
	l2 = halfedge_length[h->index],
	l3 = halfedge_length[h->next()->index];

      int idx[3];
      idx[0]=h->vertex()->index;
      idx[1]=h->next()->vertex()->index;
      idx[2]=h->prev()->vertex()->index;

      for (PetscInt j=0; j < 3; j++) idx[j] = FEM_index_map[idx[j]];

      PetscScalar mass[9], stiff[9];
      for (PetscInt j = 0; j < 9; ++j) {
	mass[j]= facet_area[i]*J4[j]/12.;
	stiff[j]= (l1*l1*J1[j]+l2*l2*J2[j]+l3*l3*J3[j])/(8.* facet_area[i]);
      }

      MatSetValues(stiff_mat, 3, idx, 3, idx, stiff, ADD_VALUES);
      MatSetValues(mass_mat, 3, idx, 3, idx, mass, ADD_VALUES);	       
      
    }

    MatAssemblyBegin(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mass_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mass_mat, MAT_FINAL_ASSEMBLY);


    clock_end();

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
    VecDestroy(&tmp);
    return inner_prod;
  }

  void TriMesh::PETSc_load_LBeigen(std::string name) {
    name.append(".ev");
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
    VecGetSize(eig_vector[0], &mat_size);
    

    // mapping FEM indices
    if (mat_size < vertex_num) {
      map_linearFEM_indices(true); 
      FEM_degree = 1;
    } else if (mat_size == vertex_num) {
      map_linearFEM_indices(false);
      FEM_degree = 1;
    } else if (mat_size < vertex_num + halfedge_num + facet_num) {
      map_cubicFEM_indices(true);
      FEM_degree = 3;
    } else if (mat_size <= vertex_num + halfedge_num + facet_num) {
      map_cubicFEM_indices(false);
      FEM_degree = 3;
    }



    // to use feature selection
    /*
    PetscScalar *vtmp_array;
    for (int i=0 ; i< fbase_size; ++i) {
      VecCopy(eig_vector[i], vtmp);
      VecGetArray(vtmp, &vtmp_array);
      for (int j = 0; j< mat_size; ++j)
	vtmp_array[j] = std::exp(- total_area * vtmp_array[j] * vtmp_array[j]/ eig_vector_sqr_norm[i] / 2.);
      VecRestoreArray(vtmp, &vtmp_array);
      std::cout<< vec_inner_prod(vtmp, vtmp) <<std::endl;
    }

    */


    clock_end();

  }

  void TriMesh::PETSc_load_vertex_eig_vector(int i, ScalarFunction & f) {
    PetscScalar *array;
    VecGetArray(eig_vector[i], &array);
    for (int j=0; j< vertex_num; ++j)
      if (FEM_index_map[j]>=0)
	f[j] = array[FEM_index_map[j]];
      else
	f[j] = 0;
  }


  void TriMesh::PETSc_export_Fourier_base(std::string name) {
    clock_start("Export Fourier base");

    std::string fbasedata_name = name; fbasedata_name.append(".fbase");
    std::ofstream fbasedata(fbasedata_name.c_str(), std::ios::out|std::ios::binary);
    Vec vtmp, vtmp2;          
    VecDuplicate(eig_vector[0], &vtmp);
    VecDuplicate(eig_vector[0], &vtmp2);

    MESHTK_SCALAR_TYPE *v;
    v = new MESHTK_SCALAR_TYPE[eig_num-1];

    for (int i= 1; i < eig_num; ++i) {
      
      VecPointwiseMult(vtmp, eig_vector[i], eig_vector[i]);      
      MatMult(mass_mat, vtmp, vtmp2);
      
      for (int j= 1; j < eig_num; ++j) {
	VecDot(vtmp2, eig_vector[j], &v[j-1]);
	v[j-1] *= std::sqrt(total_area) /eig_vector_sqr_norm[i]/std::sqrt(eig_vector_sqr_norm[j]);
      }      

      fbasedata.write((char*) v, sizeof(MESHTK_SCALAR_TYPE) * (eig_num-1));	

    }
    delete [] v;
    VecDestroy(&vtmp);
    VecDestroy(&vtmp2);

    fbasedata.close();

    clock_end();
  }






  double TriMesh::update_vertex_biharmonic_distance(int source_vertex_index,
						    ScalarFunction & biharmonic_distance) {    

    int vertex_size = biharmonic_distance.size();
    //assert(vertex_size == vertex_num);    
   
    
    for (int j=0; j < vertex_size; ++j) biharmonic_distance[j] =0;

    for (int i=eig_num-1; i> 0; --i) {
    //for (int i = 99; i>0; --i) {
      PetscScalar *vector;
      //PetscScalar *htrace;
      VecGetArray(eig_vector[i], &vector);
      //VecGetArray(bihd_trace, &htrace);

      // weighted by heat trace, probably not a good idea
      //PetscScalar target = vector[source_vertex_index] / std::sqrt(htrace[source_vertex_index]), 
      PetscScalar target =0;
      if (FEM_index_map[source_vertex_index] >=0)
	target = vector[FEM_index_map[source_vertex_index]];
      
      PetscScalar sqr_eigenvalue = eig_value[i] * eig_value[i]; 

      for (int j=0; j < vertex_size; ++j) {

	//double tmp = vector[j] / std::sqrt(htrace[j]) -target;
	double tmp = -target;
	if (FEM_index_map[j] >=0) 
	  tmp= vector[FEM_index_map[j]] - target;

	biharmonic_distance[j] += tmp * tmp / eig_vector_sqr_norm[i] / sqr_eigenvalue;
      }
    }

    double total_distance=0;
    for (int j=0; j < vertex_size; ++j) {
      // comment out for kernel matrix assembly
      //biharmonic_distance[j] = std::sqrt(biharmonic_distance[j]);
      total_distance += biharmonic_distance[j];
    }

    
    return total_distance / vertex_size;
  }

  
  int TriMesh::PETSc_assemble_export_BiHDM(std::vector<int> & sampling, //init sampling provided, keypoint based
					   int addition_size, // expect additional sampling size						   
					   std::string name) {    
    clock_start("Assemble Nystrom sampling of BiHDM");

    name += ".bihdmat";

    int total_size = assemble_export_Nystrom_matrix(sampling,
						    addition_size,
						    name,
						    &update_vertex_biharmonic_distance);
    clock_end();


    return total_size;
  }

  void TriMesh::PETSc_assemble_Fourier_BiHDM(int fbase_size) {
    clock_start("Assemble Fourier phase of BiHDM");
    PetscInt *nnz;

    if (fbase_size == 0) fbase_size = eig_num;

    // to compute hihd_trace
    Vec vtmp;
    
    VecDuplicate(eig_vector[0], &bihd_trace);
    VecSet(bihd_trace, 0.);
    VecDuplicate(eig_vector[0], &vtmp);
    
    for (int i=fbase_size-1;i>0;--i) {
      VecPointwiseMult(vtmp, eig_vector[i], eig_vector[i]);
      VecAXPY(bihd_trace, 1./(eig_value[i] * eig_value[i])/eig_vector_sqr_norm[i], vtmp);
    }
    VecDestroy(&vtmp);
    
    
    // init sparse symmtric fbihd_mat
    nnz = new PetscInt[fbase_size];
    nnz[0]=fbase_size; for (int i=1;i<fbase_size;i++) nnz[i]=1;
    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, fbase_size, fbase_size, 0, nnz, &fbihd_mat);
    delete [] nnz;

    // assembly
    //nnz = new PetscInt[fbase_size];
    //value = new PetscScalar[fbase_size];        
    PetscScalar dvalue =0., value;
    PetscInt j=0;    
    for (int i=fbase_size-1;i>0;--i){
      value=vec_inner_prod(eig_vector[i], bihd_trace) * std::sqrt(total_area)/std::sqrt(eig_vector_sqr_norm[i]);
      dvalue += 2/(eig_value[i] * eig_value[i]);
      MatSetValues(fbihd_mat, 1, &j, 1, &i, &value, INSERT_VALUES);
    }
    
    //delete [] nnz; 
    //delete [] value;
    MatSetValues(fbihd_mat, 1, &j, 1, &j, &dvalue, INSERT_VALUES);

    for (int i=fbase_size-1;i>0;--i) {
      dvalue=-2./(eig_value[i] * eig_value[i]);
      MatSetValues(fbihd_mat, 1, &i, 1, &i, &dvalue, INSERT_VALUES);
    }


    MatAssemblyBegin(fbihd_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fbihd_mat, MAT_FINAL_ASSEMBLY);

    
    clock_end();
  }

  void TriMesh::PETSc_export_Fourier_BiHDM( std::string name) {

    name += ".fbihdmat";
    PETSc_export_mat(name.c_str(), fbihd_mat);

  }



}




