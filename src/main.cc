/*
  FILE: main.cc This file is part of MeshTK.
  It is a C++ source file which implements some useful file IO by calling
  functions of meshtk library.
  
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
    TCLAP::CmdLine cmd("3D Mesh Toolkit, Supporting mesh (pre-)processing, visualization, ... \nContributor: Jianbo YE<yelpoo@gmail.com>", ' ', "0.1");

    OptString inputMeshName("i", "input_mesh_name", "File name of input mesh without file extension", true, "", "string", cmd); 
    OptString outputMeshName("", "output_mesh_name", "File name of output mesh without file extension", false, "out", "string", cmd);
    OptString inputMeshType("", "input_mesh_type", "Input mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    OptString outputMeshType("", "output_mesh_type", "output mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    
    OptInt smoothMeshIteration("s", "smooth_mesh_iter", "Number of Guassian smoothing iterations", false, 0, "unsigned int", cmd);


    OptScalar smoothMeshCoefficient("", "smooth_mesh_coeff", "Coefficient specified for Guassian smoothing", false, 1., "float", cmd);

    OptBool outputMeshSwitch("o", "output_mesh_enable", "Enable mesh Output before program exits", cmd, false) ;
    OptBool viewMeshOnly("v", "view_mesh_only", "View mesh without other rending", cmd, false);
    OptBool viewMeshCurvature("", "view_mesh_curv", "View mesh with color ramping of mean curvature", cmd, false);

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
    // mesh Gaussian smooth iteration
    for (int i=0; i<smoothMeshIteration.getValue(); ++i) {
      std::cout<< "Smooth Iteration Count: " << i << " with coefficient "<< smoothMeshCoefficient.getValue() << std::endl;
      mesh.gaussian_smooth(smoothMeshCoefficient.getValue());
      mesh.update_base();
    }


    /***************************************************************************/
    // region to test
    /*
    mesh.update_vertex_salient(5,1);
    meshtk::MeshViewer viewer(argc, argv);
    meshtk::BooleanFunction *salient_points = (meshtk::BooleanFunction *) mesh.attribute_extract(MESHTK_VERTEX_SALIENT);
    meshtk::MeshMarker painter(&mesh, salient_points);
    viewer.add_painter(&painter);
      
    viewer.init();// call this func last before loop
    viewer.view();
    */

    //std::cout<<mesh.update_vertex_neighbor(3.)<<std::endl;

    /***************************************************************************/    
    // output and display
    if (outputMeshSwitch.getValue())
      mesh.write(outputMeshName.getValue(), outputMeshType.getValue());

    if (viewMeshOnly.getValue()) {
      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshPainter painter(&mesh);
      viewer.add_painter(&painter);

      viewer.init();// call this func last before loop
      viewer.view();
    }
    else if (viewMeshCurvature.getValue()) {
      //Compute facet curvature tensor
      mesh.update_curvature();


      meshtk::MeshViewer viewer(argc, argv);
      meshtk::ScalarFunction *mean_curv = (meshtk::ScalarFunction *) mesh.attribute_extract(MESHTK_VERTEX_HCURV);
      meshtk::MeshRamper painter(&mesh, mean_curv);
      viewer.add_painter(&painter);

      viewer.init();// call this func last before loop
      viewer.view();
    }

  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







