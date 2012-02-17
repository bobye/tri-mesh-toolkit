/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2010, Universidad Politecnica de Valencia, Spain

   This file is part of SLEPc.
      
   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY 
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS 
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for 
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a generalized eigensystem Ax=kBx with matrices loaded from a file.\n"
  "This example works for both real and complex numbers.\n\n"
  "The command line options are:\n"
  "  -matrix <filename>, where <filename> = matrix (A) file in PETSc binary form.\n"
  "  -mass_matrix <filename>, where <filename> = matrix (B) file in PETSc binary form.\n"
  "  -eps_nev <num>, where <num> = number of eigenvalues\n."
  "  -shift_val <float>, <float> = shift of eigenvalues\n"
  "  -output_file <dirname>, where <dirname> = directory to store eigenpairs\n"
  "  -dirichlet, whether there exists an null space x=1.";

#include "slepceps.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char **argv )
{
  Mat         	 A,B;		  /* matrices */
  EPS         	 eps;		  /* eigenproblem solver context */
  const EPSType  type;
  ST          	 st;		  /* spectral transformation context */
  Vec            eigenvector;
  Vec            x;
  //  KSP            ksp;
  //  MatNullSpace   musp;

  PetscReal   	 error, tol, re, im;
  PetscScalar 	 kr, ki, shift_val=0.;
  PetscErrorCode ierr;
  PetscInt    	 nev, maxit, i, its, lits, nconv;
  char        	 filename[256], output_filename[256];
  PetscViewer 	 viewer,fd,fde;
  PetscBool  	 flg, output_filename_present, flgA;
  PetscBool      dirichlet;
  char                       tmp_str[PETSC_MAX_PATH_LEN];




  SlepcInitialize(&argc,&argv,(char*)0,help);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-dirichlet",&dirichlet,PETSC_NULL);CHKERRQ(ierr);


  ierr = PetscOptionsGetScalar(PETSC_NULL,"-shift_val",&shift_val,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-output_file",output_filename,256, &output_filename_present);




  ierr = PetscOptionsGetString(PETSC_NULL,"-mass_matrix",filename,256,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem stored in file.\n\n");CHKERRQ(ierr);

    //    SETERRQ(1,"Must indicate a file name for matrix B with the -mass_matrix option.");
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    ierr = MatLoad(B, viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetString(PETSC_NULL,"-matrix",filename,256,&flgA);CHKERRQ(ierr);
  if (!flgA) {
    SETERRQ(PETSC_COMM_WORLD, 1,"Must indicate a file name for matrix A with the -matrix option.");
  }

#if defined(PETSC_USE_COMPLEX)
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...");CHKERRQ(ierr);
#else
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...");CHKERRQ(ierr);
#endif
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  if (flgA) {
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr); 
    ierr = MatLoad(A, viewer);CHKERRQ(ierr);}

  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  
  

  // ierr = PetscOptionsGetInt(PETSC_NULL,"-n_eigs",&n_eigs,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t[done]\n");CHKERRQ(ierr);  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


  /* 
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /* 
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  if (flg) {
    ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRQ(ierr);
  }
  else {
    ierr = EPSSetOperators(eps,A,PETSC_NULL);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  }
  

 
  /* 

  */


  //ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);
  
  //ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  if (flg) {
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = EPSSetTarget(eps,-shift_val);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
    ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  }
  else{
    //ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE); CHKERRQ(ierr);
  }



  ierr = MatGetVecs(A,&x,PETSC_NULL);CHKERRQ(ierr);
  if (flg){
    // COMMENT: for Dirichlet boundary problem
    if (!dirichlet) {
      ierr = VecSet(x,1.0);CHKERRQ(ierr);
      ierr = EPSSetDeflationSpace(eps,1,&x);CHKERRQ(ierr);
    }

  }


  //  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 1, &x, &musp); CHKERRQ(ierr);
  
  //  ierr = STGetKSP(st, &ksp);CHKERRQ(ierr);
  //  ierr = KSPSetNullSpace(ksp, musp);CHKERRQ(ierr);
  //  ierr = STSetKSP(st, ksp);CHKERRQ(ierr);



  //
  







  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  //  ierr = EPSGetDimensions(eps,&n_eigs, PETSC_NULL ,PETSC_NULL);CHKERRQ(ierr);



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
  ierr = EPSGetOperationCounters(eps,PETSC_NULL,PETSC_NULL,&lits);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Get number of converged eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);

  sprintf( tmp_str, "%s/_ev.ascii", output_filename);
  if(output_filename_present)   PetscViewerASCIIOpen(PETSC_COMM_WORLD, tmp_str, &fde);

  if (output_filename_present&&flg&& !dirichlet)
    {
      
      ierr = PetscViewerASCIIPrintf(fde,"%e\n",0);	   
      
       sprintf( tmp_str, "%s/_%d.petsc", output_filename, 0 );
       PetscViewerBinaryOpen(PETSC_COMM_WORLD, tmp_str, FILE_MODE_WRITE, &fd);
	  /* PetscViewerSetFormat(fd,PETSC_VIEWER_ASCII_MATLAB); */
       ierr = VecView(x,fd); CHKERRQ(ierr);
       ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    }



  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k             ||Ax-kBx||/||kx||\n"
         "  --------------------- ------------------\n" );CHKERRQ(ierr);
    for( i=0; i<nconv; i++ ) {
      /* 
         Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
         ki (imaginary part)
      */
      VecDuplicate(x, &eigenvector);

      ierr = EPSGetEigenpair(eps,i,&kr,&ki,eigenvector,PETSC_NULL);CHKERRQ(ierr);

      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      //      eigs[i]=re;

      if( im != 0.0 ) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," % 6g %+6g i",re,im);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  % 12.12e      ",re); CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);



      if (output_filename_present)
	{
	  
	  ierr = PetscViewerASCIIPrintf(fde,"%e\n",re);	   
	  
	  if (dirichlet) sprintf( tmp_str, "%s/_%d.petsc", output_filename, i );
	  else sprintf( tmp_str, "%s/_%d.petsc", output_filename, i+1 );

	  PetscViewerBinaryOpen(PETSC_COMM_WORLD, tmp_str, FILE_MODE_WRITE, &fd);
	  /* PetscViewerSetFormat(fd,PETSC_VIEWER_ASCII_MATLAB); */
	  ierr = VecView(eigenvector,fd); CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	}
            
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
  }

  if(output_filename_present) ierr = PetscViewerDestroy(&fde);CHKERRQ(ierr);






  
  /* 
     Free work space
  */
  
  //   VecDestroyVecs(eigenvectors, n_eigs);
   //free(eigs);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&x); CHKERRQ(ierr); 
  ierr = VecDestroy(&eigenvector); CHKERRQ(ierr); 
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  if (flg) {ierr = MatDestroy(&B); CHKERRQ(ierr);}
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;
}



