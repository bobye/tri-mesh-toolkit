/*

  Lapack_wrapper.h
  Header file for Lapack wrapper functions. See lapack_wrapper.c for more details
  Rob Heylen, aug 2006
  http://itf.fys.kuleuven.be/~rob/computer/lapack_wrapper/index.htm

*/

#ifndef _LAPACK_WRAPPER_H_
#define _LAPACK_WRAPPER_H_

//namespace lapack_wrapper {
#ifdef __cplusplus
extern "C" {
#endif


extern void mat2cvec(int m, int n, double** mat, double* cvec);

extern void mat2fvec(int m, int n, double** mat, double* fvec);

extern void cvec2mat(int m, int n, double* cvec, double** mat);

extern void fvec2mat(int m, int n, double* fvec, double** mat);

extern void cvec2fvec(int m, int n, double *cvec, double* fvec);

extern void fvec2cvec(int m, int n, double *fvec, double* cvec);

extern void matrix_matrix_mult(int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);

extern void matrix_add(int m, int n, double a, double **X, double **Y);

extern void vector_add(int n, double a, double *X, double *Y);

extern double dotprod(int n, double *X, double *Y);

extern void vector_copy(int n, double* X, double* Y);

extern void vector_scale(int n, double a, double* X);

extern void matrix_vector_mult(int m, int n, double a, double b, double *A, double *x, double *y);

extern void matrix_transpose(int m, int n, double *X, double *Y);

extern int eigen_decomposition(int n, double* X, double *eigvec, double *eigval);

extern int matrix_square_root_n(int n, double *X, double *I, double *Y);

extern int matrix_square_root(int n, double *X, double *Y);

extern int matrix_invert(int n, double *X, double *Y);

extern int linsolve(int n, double *A, double *B, double *x);

//}
#ifdef __cplusplus
}
#endif

#endif /* _LAPACK_WRAPPER_H_ */
