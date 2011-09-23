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
typedef TCLAP::ValueArg<std::string>                 Opt_String;
typedef TCLAP::ValueArg<double>                      Opt_scalar;
typedef TCLAP::ValueArg<int>                         Opt_Int;
typedef TCLAP::SwitchArg                             Opt_Bool;

//#include "mesh_topo.hh"
#include "meshtk/TriMesh.hh"
#include "meshtk/MeshViewer.hh"



int main(int argc, char** argv){

  try {
    TCLAP::CmdLine cmd("3D Mesh Toolkit, Supporting mesh (pre-)processing, visualization, ... \nContributor: Jianbo YE<yelpoo@gmail.com>", ' ', "0.1");

    Opt_String inputMeshName("i", "input_mesh", "File name of input mesh without file extension", true, "", "string", cmd); 
    Opt_String outputMeshName("o", "output_mesh", "File name of output mesh without file extension", false, "out", "string", cmd);
    Opt_String inputMeshType("", "input_mesh_type", "Input mesh file format, candidates are off(default), ...", false, "off", "string", cmd);
    Opt_String outputMeshType("", "output_mesh_type", "output mesh file format, candidates are off(default), ...", false, "off", "string", cmd);



    cmd.parse( argc, argv );

    std::string input_mesh_name  = inputMeshName.getValue();
    std::string output_mesh_name = outputMeshName.getValue();
    std::string input_mesh_type  = inputMeshType.getValue();
    std::string output_mesh_type = outputMeshType.getValue();

    
    meshtk::TriMesh mesh;// index mapping from array to halfedge data structure

    mesh.read(input_mesh_name, input_mesh_type);

    mesh.init_index();
    
    // Update mesh infomation
    mesh.update_base();
    

    //mesh processing
    /***************************************************************************/

    //Compute facet curvature tensor
    mesh.update_curvature();
    /***************************************************************************/    


    //mesh.write(output_mesh_name, output_mesh_type);
    
    meshtk::MeshViewer viewer(argc, argv);
    //meshtk::MeshPainter painter(&mesh);
    meshtk::ScalarFunction *mean_curv = (meshtk::ScalarFunction *) mesh.attribute_extract(MESHTK_VERTEX_HCURV);
    meshtk::MeshRamper painter(&mesh, mean_curv);

    viewer.add_painter(&painter);

    viewer.init();// call this func last before loop
    viewer.view();

  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







