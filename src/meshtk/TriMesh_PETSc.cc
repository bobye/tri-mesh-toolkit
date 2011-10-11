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

  Mat mass_mat, stiff_mat;


  void TriMesh::PETSc_assemble_cubicFEM_LBmat(){
    clock_start("Cubic elements assembly");

    // matrix initialization
    PetscInt n = vertex_num + halfedge_num + facet_num;

    PetscInt *nnz=new PetscInt[n];

    for (PetscInt i=0;i<vertex_num;i++)
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


    clock_end();
  }

  void TriMesh::PETSc_assemble_linearFEM_LBmat(){
  }




}

