/*
  FILE: meshtk/TriMesh_PETSc.cc This file is part of MeshTK.
  It is a C++ header file which implement PETSc matrix for mesh.
  
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

#include "meshtk/TriMesh.hh"
#include "meshtk/mesh_assist.hh"
#include <petscmat.h>


namespace meshtk {

  static char help[] = "";

  void TriMesh::PETSc_init(int argc, char **argv){
    PetscInitialize(&argc,&argv,(char *)0,help);
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



  void TriMesh::PETSc_assemble_cubicFEM_LBmat(std::string name){
    clock_start("Cubic elements assembly");
    Mat mass_mat, stiff_mat;

    // matrix initialization
    PetscInt n = vertex_num + halfedge_num + facet_num;

    PetscInt *nnz=new PetscInt[n];

    for (PetscInt i=0;i<vertex_num;++i)
      nnz[i]=IV[i]->vertex_degree()*6+1;

    for (PetscInt i=0;i<halfedge_num;i++)
      nnz[i+vertex_num]=12;

    for (PetscInt i=0;i<facet_num;i++)
      nnz[i+vertex_num+halfedge_num]=1;

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, n, n, 0, nnz, &mass_mat);
    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, n, n, 0, nnz, &stiff_mat);

    MatSetOption(mass_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
    MatSetOption(stiff_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    delete nnz;

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

    PetscViewer fmass, fstiff;
    std::string filename_stiff = name, filename_mass = name;
    filename_stiff += ".stiff"; filename_mass += ".mass";

    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename_mass.c_str() , FILE_MODE_WRITE, &fmass);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename_stiff.c_str(), FILE_MODE_WRITE, &fstiff);
      
    MatView(mass_mat, fmass);
    MatView(stiff_mat, fstiff);

    PetscViewerDestroy(&fmass); PetscViewerDestroy(&fstiff);
    MatDestroy(&mass_mat); MatDestroy(&stiff_mat);


    clock_end();
  }

  void TriMesh::PETSc_assemble_linearFEM_LBmat( std::string name){


  }




}

