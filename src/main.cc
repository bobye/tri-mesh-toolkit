/*
  FILE: MeshTK.cc This file is part of tri-mesh-toolkit.
  It is a C++ source file which implements some useful file IO by calling
  functions of meshtk library.
  
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


#include <string>
#include <algorithm>

#include <tclap/CmdLine.h>
typedef TCLAP::ValueArg<std::string>                 OptString;
typedef TCLAP::ValueArg<double>                      OptScalar;
typedef TCLAP::ValueArg<int>                         OptInt;
typedef TCLAP::SwitchArg                             OptBool;

#include <meshtk/TriMesh.hh>
#include <meshtk/DynamicTriMesh.hh>
#include <meshtk/MeshViewer.hh>



int main(int argc, char** argv){

  try {
    TCLAP::CmdLine cmd("3D Mesh Toolkit, supporting mesh (pre-)processing, shape analysis, visualization, ... \nContributor: Jianbo YE<yelpoo@gmail.com> \nhttp://code.google.com/p/tri-mesh-toolkit/", ' ', MESHTK_VERSION);




    OptInt smoothMeshIteration("s", "smooth_mesh_iter", "Number of Guassian smoothing iterations", false, 0, "unsigned int", cmd);
    OptInt viewMeshGeodesicDist("", "view_geodesic_source", "View geodesic distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshBiharmonicDist("", "view_biharmonic_source", "View biharmonic distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshEigenvector("", "view_eigenvector", "View color ramping of the i-th eigenvector", false, -1, "unsigned int", cmd);
    OptInt removeMeshFacetInt("", "remove_facets_int","Remove n facets randomly", false, -1, "unsigned int", cmd);




    OptScalar smoothMeshCoefficient("", "smooth_mesh_coeff", "Coefficient specified for Guassian smoothing", false, 1., "float", cmd);
    OptScalar addMeshNoise("", "noise_mesh", "Coefficient specified for uniform mesh noise added", false, 0., "float", cmd);
    OptScalar removeMeshFacetPTG("", "remove_facets_ptg", "Remove p%% facets randomly", false, 0., "float", cmd);

    

    OptBool exportPCAviewGCL("", "export_PCA_view_gcl", "Export Geomview Command Language to view model by PCA analysis", cmd, false);

    OptString loadMeshCurvName("", "load_mesh_curv_name", "load vertex curvature from file with specified name prefix", false, "", "string", cmd);
    OptBool loadMeshCurv("", "load_mesh_curv", "load vertex curvature from file", cmd, false);



    OptString loadMeshLBeigenName("", "load_FEM_LBeigen_name", "load FEM eigenvalues and eigenvectors of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptBool loadMeshLBeigen("g", "load_FEM_LBeigen", "load FEM eigenvalues and eigenvectors of Laplace Beltrami operator", cmd, false);

    OptString loadMeshLBmatName("", "load_FEM_LBmat_name", "Load FEM mass and stiff matrix of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptBool loadMeshLBmat("m", "load_FEM_LBmat", "Load FEM mass and stiff matrix of Laplace Beltrami operator", cmd, false);


    OptString exportMeshFBiHDMName("", "export_Fourier_BiHDM_name", "Export Fourier phase of square Biharmonic distance matrix with specified name prefix", false, "", "string", cmd);
    OptInt exportMeshFBiHDMSize("", "export_Fourier_BiHDM_size", "Export Fourier phase of square Biharmonic distance matrix with specified size of bases", false, 0, "unsigned int", cmd);
    OptBool exportMeshFBiHDM("f", "export_Fourier_BiHDM", "Export Fourier phase of square Biharmonic distance matrix", cmd, false);


    OptString exportMeshFBaseName("", "export_Fourier_base_name", "Export Fourier base square matrix with specified name prefix", false, "", "string", cmd);
    OptBool exportMeshFBase("b", "export_Fourier_base", "Export Fourier base square matrix", cmd, false);

    OptBool enableDirichletBC("", "Dirichlet", "Enable Dirichlet boundary condition for FEM formulation", cmd, false);


    OptString exportMeshCubicLBmatName("", "export_cubic_FEM_LBmat_name", "Export cubic FEM mass and stiff matrix of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptBool exportMeshCubicLBmat("c", "export_cubic_FEM_LBmat", "Export cubic FEM mass and stiff matrix of Laplace Beltrami operator", cmd, false);

    OptString exportMeshLinearLBmatName("", "export_linear_FEM_LBmat_name", "Export linear FEM mass and stiff matrix of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptBool exportMeshLinearLBmat("l", "export_linear_FEM_LBmat", "Export linear FEM mass and stiff matrix of Laplace Beltrami operator", cmd, false);



    OptBool viewMeshOnly("v", "view_mesh_only", "View mesh without other rending", cmd, false);
    OptBool viewMeshCurvature("", "view_mesh_curv", "View mesh with color ramping of mean curvature", cmd, false);


    OptString outputMeshType("", "output_mesh_type", "output mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    OptString outputMeshName("", "output_mesh_name", "File name of output mesh without file extension", false, "out", "string", cmd);
    OptBool outputMeshSwitch("o", "output_mesh_enable", "Enable mesh Output before program exits", cmd, false) ;


    OptString inputMeshType("", "input_mesh_type", "Input mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    OptString inputMeshName("i", "input_mesh_name", "File name of input mesh without file extension", true, "", "string", cmd); 

    // process input argument
    cmd.parse( argc, argv );

    
    /////////////////////////////////////////////////////////////////////////////////    
    meshtk::DynamicTriMesh mesh;// index mapping from array to halfedge data structure
    // load mesh by name and type, required for program run successfully
    mesh.read(inputMeshName.getValue(), inputMeshType.getValue());
    mesh.init_index();// initialize mesh    
    // Update mesh infomation, prompt after init_index()
    mesh.update_base();
    /////////////////////////////////////////////////////////////////////////////////
    // dynamic mesh processing starts here
    // add mesh noise
    if (addMeshNoise.getValue() > 0) {
      mesh.add_mesh_noise(addMeshNoise.getValue());
      mesh.update_base();
    }

    // remove facets randomly
    if (removeMeshFacetInt.getValue() > 0 && removeMeshFacetPTG.getValue() == 0) {
      mesh.remove_mesh_facets(removeMeshFacetInt.getValue());
      mesh.init_index();
      mesh.update_base();
    } else if (removeMeshFacetInt.getValue() > 0) {
      mesh.remove_mesh_facets(removeMeshFacetInt.getValue(), removeMeshFacetPTG.getValue());
      mesh.init_index();
      mesh.update_base();
    } else {
      mesh.remove_mesh_facets(-1, removeMeshFacetPTG.getValue());
      mesh.init_index();
      mesh.update_base();
    }

    // mesh Gaussian smooth iteration
    // not able to handle border now
    for (int i=0; i<smoothMeshIteration.getValue(); ++i) {
      std::cout<< "Smooth Iteration Count: " << i << " with coefficient "<< smoothMeshCoefficient.getValue() << std::endl;
      
      mesh.gaussian_smooth(smoothMeshCoefficient.getValue());
      mesh.update_base();
    }




    if (exportPCAviewGCL.getValue())
      mesh.PCAview2gcl(inputMeshName.getValue());
    


    /////////////////////////////////////////////////////////////////////////////////
    // shape analysis starts here
    mesh.PETSc_init(argc, argv);

    if (exportMeshLinearLBmat.getValue()) {
      mesh.PETSc_assemble_linearFEM_LBmat(enableDirichletBC.getValue());

      if (exportMeshLinearLBmatName.getValue().compare("") == 0) 
	mesh.PETSc_export_LBmat(inputMeshName.getValue());
      else mesh.PETSc_export_LBmat(exportMeshLinearLBmatName.getValue());
      
    }
    else if (exportMeshCubicLBmat.getValue()) {
      mesh.PETSc_assemble_cubicFEM_LBmat(enableDirichletBC.getValue());

      if (exportMeshCubicLBmatName.getValue().compare("") == 0) 
	mesh.PETSc_export_LBmat(inputMeshName.getValue());
      else mesh.PETSc_export_LBmat(exportMeshCubicLBmatName.getValue());
      
    }
    else if (loadMeshLBmat.getValue()) {
      if (loadMeshLBmatName.getValue().compare("") == 0) 
	mesh.PETSc_load_LBmat(inputMeshName.getValue());
      else mesh.PETSc_load_LBmat(loadMeshLBmatName.getValue());
    }
    
    if (loadMeshLBeigen.getValue()) {
      if (loadMeshLBeigenName.getValue().compare("") == 0)
	mesh.PETSc_load_LBeigen(inputMeshName.getValue());
      else mesh.PETSc_load_LBeigen(loadMeshLBeigenName.getValue());      
    }

    if (exportMeshFBiHDM.getValue() && loadMeshLBeigen.getValue()) {
      mesh.PETSc_assemble_Fourier_BiHDM(exportMeshFBiHDMSize.getValue());
      if (exportMeshFBiHDMName.getValue().compare("") == 0)
	mesh.PETSc_export_Fourier_BiHDM(inputMeshName.getValue());
      else 
	mesh.PETSc_export_Fourier_BiHDM(exportMeshFBiHDMName.getValue());
    }


    if (exportMeshFBase.getValue() && loadMeshLBeigen.getValue()) {
      if (exportMeshFBaseName.getValue().compare("") == 0)
	mesh.PETSc_export_Fourier_base(inputMeshName.getValue());
      else mesh.PETSc_export_Fourier_base(exportMeshFBaseName.getValue());
    }








    /***************************************************************************/    
    // region to test
    //mesh.detect_vertex_salient(10,1);

    /***************************************************************************/    
    // output and display

    if (outputMeshSwitch.getValue())
      mesh.write(outputMeshName.getValue(), outputMeshType.getValue());

    if (viewMeshOnly.getValue()) {
      mesh.update_compact_base();
      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshPainter painter(&mesh);
      viewer.add_painter(&painter);

      viewer.init();// call this func last before loop
      viewer.view();
    }
    else if (viewMeshCurvature.getValue()) {
      //Compute facet curvature tensor
      mesh.update_compact_base();
      mesh.update_curvature();

      meshtk::MeshViewer viewer(argc, argv);
      meshtk::ScalarFunction *mean_curv = (meshtk::ScalarFunction *) mesh.attribute_extract(MESHTK_VERTEX_HCURV);
      meshtk::MeshRamper painter(&mesh, mean_curv, MESHTK_COLOR_HIST);
      viewer.add_painter(&painter);

      viewer.init();// call this func last before loop
      viewer.view();
    }
    else if (viewMeshGeodesicDist.getValue()>=0) {
      unsigned GEODESIC_DIST_REG = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
      meshtk::ScalarFunction *geodesic_distance = (meshtk::ScalarFunction *) mesh.attribute_extract(GEODESIC_DIST_REG);
      mesh.update_compact_base();
      mesh.update_vertex_geodesic_distance(viewMeshGeodesicDist.getValue(), *geodesic_distance);
      
      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshRamper painter(&mesh, geodesic_distance, MESHTK_COLOR_CONT);
      viewer.add_painter(&painter);

      viewer.init();
      viewer.view();      
    }
    else if (viewMeshBiharmonicDist.getValue() >=0 && loadMeshLBeigen.getValue()) {
      unsigned BIHARMONIC_DIST_REG = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
      meshtk::ScalarFunction *biharmonic_distance = (meshtk::ScalarFunction *) mesh.attribute_extract(BIHARMONIC_DIST_REG);
      mesh.update_compact_base();
      mesh.update_vertex_biharmonic_distance(viewMeshBiharmonicDist.getValue(), *biharmonic_distance);

      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshRamper painter(&mesh, biharmonic_distance, MESHTK_COLOR_CONT);
      viewer.add_painter(&painter);

      viewer.init();
      viewer.view();

    }
    else if (viewMeshEigenvector.getValue() >=0 && loadMeshLBeigen.getValue()) {
      unsigned EIGEN_VECTOR = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
      meshtk::ScalarFunction * eigenvector_scalar = (meshtk::ScalarFunction *) mesh.attribute_extract(EIGEN_VECTOR);
      mesh.update_compact_base();
      mesh.PETSc_load_vertex_eig_vector(viewMeshEigenvector.getValue(), *eigenvector_scalar);
      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshRamper painter(&mesh, eigenvector_scalar, MESHTK_COLOR_CONT);
      viewer.add_painter(&painter);

      viewer.init();
      viewer.view();
    }

    /***************************************************************************/
    // region to destroy

    mesh.PETSc_destroy();



  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}








