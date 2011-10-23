/*
  FILE: main.cc This file is part of tri-mesh-toolkit.
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

    OptString inputMeshName("i", "input_mesh_name", "File name of input mesh without file extension", true, "", "string", cmd); 
    OptString outputMeshName("", "output_mesh_name", "File name of output mesh without file extension", false, "out", "string", cmd);
    OptString inputMeshType("", "input_mesh_type", "Input mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    OptString outputMeshType("", "output_mesh_type", "output mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    OptString exportMeshCubicLBmatName("", "export_cubic_FEM_LBmat_name", "Export cubic FEM mass and stiff matrix of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptString exportMeshFBiHDMName("", "export_Fourier_BiHDM_name", "Export Fourier phase of square Biharmonic distance matrix with specified name prefix", false, "", "string", cmd);
    OptString loadMeshLBmatName("", "load_FEM_LBmat_name", "Load FEM mass and stiff matrix of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptString loadMeshLBeigenName("", "load_FEM_LBeigen_name", "load FEM eigenvalues and eigenvectors of Laplace Beltrami operator with specified name prefix", false, "", "string", cmd);
    OptString loadMeshCurvName("", "load_mesh_curv_name", "load vertex curvature from file with specified name prefix", false, "", "string", cmd);
    OptString loadMeshSIFTName("", "load_mesh_SIFT_name", "Load mesh local descriptors for all vertices with specified name prefix", false, "", "string", cmd);
    OptString exportMeshSIFTName("", "export_mesh_SIFT_name", "Export mesh local descriptors for all vertices with specified name prefix", false, "", "string", cmd);
    OptString loadMeshHKSName("", "load_mesh_HKS_name", "Load mesh HKS descriptors for all vertices with specified name prefix", false, "", "string", cmd);
    OptString exportMeshHKSName("", "export_mesh_HKS_name", "Export mesh HKS descriptors for all vertices with specified name prefix", false, "", "string", cmd);
    OptString exportMeshBiHSIFTmixDMName("", "export_BiH_SIFT_mixDM_name", "Export Nystrom sampling matrix block of biharmonic-SIFT mixed distance matrix with specified name prefix", false, "", "string", cmd);
    

    OptInt smoothMeshIteration("s", "smooth_mesh_iter", "Number of Guassian smoothing iterations", false, 0, "unsigned int", cmd);
    OptInt viewMeshGeodesicDist("", "view_geodesic_source", "View geodesic distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshBiharmonicDist("", "view_biharmonic_source", "View biharmonic distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshSIFTDist("", "view_SIFT_source", "View SIFT distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshHKSDist("", "view_HKS_source", "View HKS distance from a source vertex on mesh", false, -1, "index", cmd);
    OptInt viewMeshEigenvector("", "view_eigenvector", "View color ramping of the i-th eigenvector", false, -1, "unsigned int", cmd);




    OptScalar smoothMeshCoefficient("", "smooth_mesh_coeff", "Coefficient specified for Guassian smoothing", false, 1., "float", cmd);
    OptScalar addMeshNoise("", "noise_mesh", "Coefficient specified for uniform mesh noise added", false, 0., "float", cmd);


    
    OptBool outputMeshSwitch("o", "output_mesh_enable", "Enable mesh Output before program exits", cmd, false) ;
    OptBool viewMeshOnly("v", "view_mesh_only", "View mesh without other rending", cmd, false);
    OptBool viewMeshCurvature("", "view_mesh_curv", "View mesh with color ramping of mean curvature", cmd, false);
    OptBool exportMeshCubicLBmat("l", "export_cubic_FEM_LBmat", "Export cubic FEM mass and stiff matrix of Laplace Beltrami operator", cmd, false);
    OptBool exportMeshFBiHDM("f", "export_Fourier_BiHDM", "Export Fourier phase of square Biharmonic distance matrix", cmd, false);
    OptBool loadMeshLBmat("", "load_FEM_LBmat", "Load FEM mass and stiff matrix of Laplace Beltrami operator", cmd, false);
    OptBool loadMeshLBeigen("", "load_FEM_LBeigen", "load FEM eigenvalues and eigenvectors of Laplace Beltrami operator", cmd, false);
    OptBool loadMeshCurv("", "load_mesh_curv", "load vertex curvature from file", cmd, false);
    OptBool loadMeshSIFT("", "load_mesh_SIFT", "Load mesh local descriptors for all vertices", cmd, false);
    OptBool exportMeshSIFT("d", "export_mesh_SIFT", "Export mesh local descriptors for all vertices", cmd, false);
    OptBool loadMeshHKS("", "load_mesh_HKS", "Load mesh HKS descriptors for all vertices", cmd, false);
    OptBool exportMeshHKS("k", "export_mesh_HKS", "Export mesh HKS descriptors for all vertices", cmd, false);
    OptBool exportMeshBiHSIFTmixDM("m", "export_BiH_SIFT_mixDM", "Export Nystrom sampling matrix block of biharmonic-SIFT mixed distance matrix", cmd, false);
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

    // mesh Gaussian smooth iteration
    for (int i=0; i<smoothMeshIteration.getValue(); ++i) {
      std::cout<< "Smooth Iteration Count: " << i << " with coefficient "<< smoothMeshCoefficient.getValue() << std::endl;
      
      mesh.gaussian_smooth(smoothMeshCoefficient.getValue());
      mesh.update_base();
    }

    /////////////////////////////////////////////////////////////////////////////////
    // shape analysis starts here
    mesh.PETSc_init(argc, argv);

    if (exportMeshCubicLBmat.getValue()) {
      mesh.PETSc_assemble_cubicFEM_LBmat();
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
    
    if (exportMeshSIFT.getValue()) {
      mesh.update_compact_base();

      if (loadMeshCurv.getValue()) {
	if (loadMeshCurvName.getValue().compare("") == 0)
	  mesh.load_vertex_curvature(inputMeshName.getValue());
	else 
	  mesh.load_vertex_curvature(loadMeshCurvName.getValue());
      }	
      else
	mesh.update_curvature();
      

      mesh.update_all_vertices_SIFT();
      if (exportMeshSIFTName.getValue().compare("") == 0) 
	mesh.export_keypoint_SIFT(inputMeshName.getValue());      
      else mesh.export_keypoint_SIFT(exportMeshSIFTName.getValue());
    }
    else if (loadMeshSIFT.getValue()) {
      if (loadMeshSIFTName.getValue().compare("") == 0)
	mesh.load_all_vertices_SIFT(inputMeshName.getValue());
      else mesh.load_all_vertices_SIFT(loadMeshSIFTName.getValue());
    }

    if (exportMeshHKS.getValue() && loadMeshLBeigen.getValue()) {
      if (exportMeshHKSName.getValue().compare("") ==0)
	mesh.update_export_all_vertices_HKS(inputMeshName.getValue());
      else 
	mesh.update_export_all_vertices_HKS(exportMeshHKSName.getValue());
    }
    else if (loadMeshHKS.getValue()) {
      if (loadMeshHKSName.getValue().compare("") == 0)
	mesh.load_all_vertices_HKS(inputMeshName.getValue());
      else 
	mesh.load_all_vertices_HKS(loadMeshHKSName.getValue());
    }

    if (exportMeshFBiHDM.getValue() && loadMeshLBeigen.getValue()) {
      mesh.PETSc_assemble_Fourier_BiHDM();
      if (exportMeshFBiHDMName.getValue().compare("") == 0)
	mesh.PETSc_export_Fourier_BiHDM(inputMeshName.getValue());
      else 
	mesh.PETSc_export_Fourier_BiHDM(exportMeshFBiHDMName.getValue());
    }

    if (exportMeshBiHSIFTmixDM.getValue() && loadMeshLBeigen.getValue()) {
      mesh.update_compact_base();
      unsigned USER_MESH_KEYPOINT = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_BOOLEAN);
      meshtk::BooleanFunction *mesh_keypoint = (meshtk::BooleanFunction *) mesh.attribute_extract(USER_MESH_KEYPOINT);
  
      mesh.update_curvature();
      meshtk::ScalarFunction *mesh_hcurv = (meshtk::ScalarFunction *) mesh.attribute_extract(MESHTK_VERTEX_HCURV);
      mesh.detect_vertex_keypoint(*mesh_hcurv, *mesh_keypoint, 100);
      mesh.threshold_keypoint(.618);

      std::vector<int> keypoint_threshold_index;
      mesh.export_keypoint_index(keypoint_threshold_index);

      //mesh.PETSc_load_LBmat(inputMeshName.getValue());
      //mesh.PETSc_load_LBeigen(inputMeshName.getValue());
      /*
      if (loadMeshSIFTName.getValue().compare("") == 0)
	mesh.load_all_vertices_SIFT(inputMeshName.getValue());
      else mesh.load_all_vertices_SIFT(loadMeshSIFTName.getValue());
      */

      if (exportMeshBiHSIFTmixDMName.getValue().compare("") == 0)	
	mesh.PETSc_assemble_export_BiH_SIFTmixDM(keypoint_threshold_index, 200, inputMeshName.getValue());
      else mesh.PETSc_assemble_export_BiH_SIFTmixDM(keypoint_threshold_index, 200, exportMeshBiHSIFTmixDMName.getValue());
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
      meshtk::MeshRamper painter(&mesh, mean_curv, true);
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
      meshtk::MeshRamper painter(&mesh, geodesic_distance);
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
      meshtk::MeshRamper painter(&mesh, biharmonic_distance);
      viewer.add_painter(&painter);

      viewer.init();
      viewer.view();

    }
    else if (viewMeshSIFTDist.getValue() >=0 && loadMeshSIFT.getValue()) {
      unsigned SIFT_DIST_REG = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
      meshtk::ScalarFunction *SIFT_distance = (meshtk::ScalarFunction *) mesh.attribute_extract(SIFT_DIST_REG);
      mesh.update_compact_base();
      mesh.update_vertex_SIFT_distance(viewMeshSIFTDist.getValue(), *SIFT_distance);

      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshRamper painter(&mesh, SIFT_distance);
      viewer.add_painter(&painter);

      viewer.init();
      viewer.view();

    }
    else if (viewMeshHKSDist.getValue() >=0 && loadMeshHKS.getValue()) {
      unsigned HKS_DIST_REG = mesh.attribute_allocate(MESHTK_VERTEX, MESHTK_SCALAR);
      meshtk::ScalarFunction *HKS_distance = (meshtk::ScalarFunction *) mesh.attribute_extract(HKS_DIST_REG);
      mesh.update_compact_base();
      mesh.update_vertex_HKS_distance(viewMeshHKSDist.getValue(), *HKS_distance);

      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshRamper painter(&mesh, HKS_distance);
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
      meshtk::MeshRamper painter(&mesh, eigenvector_scalar);
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







