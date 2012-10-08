#include "meshtk/ManifoldTriMesh.hh"
#include "meshtk/mesh_assist.hh"
#include <petscmat.h>
#include <iostream>
#include <fstream>

namespace meshtk {
  extern Mat stiff_mat, mass_mat;
  extern PetscInt mat_size;


  extern std::vector<double> eig_value;
  extern std::vector<Vec> eig_vector;
  extern std::vector<double> eig_vector_sqr_norm;
  // number of eigenvalues computed
  extern PetscInt eig_num;

  void ManifoldTriMesh::PETSc_assemble_graphcut() {
    clock_start("Assembly GraphCut Matrix");
    mat_size = map_linearFEM_indices(false);

    PetscInt *nnz =new PetscInt[mat_size];
    for (PetscInt i=0;i<vertex_num;++i) nnz[i]=IV[i]->vertex_degree()+1;  

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &mass_mat);
    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, mat_size, mat_size, 0, nnz, &stiff_mat);

    delete [] nnz;
    MatSetOption(mass_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
    MatSetOption(stiff_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    ScalarFunction halfedge_rotate_distance;
    double avg_edge_rotate_distance = 0;
    halfedge_rotate_distance.resize(halfedge_num);

    int m = vertex_rotate_exp.size();
    for (PetscInt i = 0; i < halfedge_num; ++i) {
      Halfedge_handle h = IH[i];
      halfedge_rotate_distance[i] = 0;
      for (int j = 0; j < m; ++j) {
	Vector v = vertex_rotate_exp[j][h->vertex()->index] - vertex_rotate_exp[j][h->prev()->vertex()->index];
	halfedge_rotate_distance[i] +=  v*v;
      }

      avg_edge_rotate_distance += (halfedge_rotate_distance[i] = CGAL::sqrt(halfedge_rotate_distance[i]/m));
    }

    avg_edge_rotate_distance /= halfedge_num;


    for (PetscInt i = 0; i < vertex_num; ++i) {
      HV_circulator hv = IV[i]->vertex_begin();
      double total = 0; int vi = hv->vertex()->index;
      do {
	// assembly mass and stiff	
	double value =  - std::exp(- halfedge_length[hv->index] * halfedge_length[hv->index] / (2 * avg_edge_len * avg_edge_len)) * std::exp( - halfedge_rotate_distance[hv->index] * halfedge_rotate_distance[hv->index] / (2 * avg_edge_rotate_distance * avg_edge_rotate_distance));

	MatSetValues(stiff_mat, 1, &(vi), 1, &(hv->prev()->vertex()->index), &value, INSERT_VALUES);
	total -= value;
      } while (++hv != IV[i]->vertex_begin()) ;

      MatSetValues(mass_mat, 1, &vi, 1, &vi, &total, INSERT_VALUES);	       	
      MatSetValues(stiff_mat, 1, &vi, 1, &vi, &total, INSERT_VALUES);	       	

    }

    MatAssemblyBegin(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mass_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mass_mat, MAT_FINAL_ASSEMBLY);

    clock_end();
  }

  void ManifoldTriMesh::PETSc_export_graphcut_vectors(std::string name){
    clock_start("Export GraphCut Vectors");
    name.append(".embd");
    std::ofstream fid(name.c_str());
    
    fid << vertex_num << " " << eig_num -1 << std::endl;

    std::vector<PetscScalar *> vector; vector.resize(eig_num);
    
    for (int j=0; j<eig_num; ++j) VecGetArray(eig_vector[j], &vector[j]);

    
    for (int i = 0; i < vertex_num; ++i) {
      for (int j = 1; j < eig_num; ++j) 
	fid << vector[j][i]/std::sqrt(eig_value[j]) << " ";
      fid << std::endl;
    }

    fid.close();
    clock_end();
  }
}
