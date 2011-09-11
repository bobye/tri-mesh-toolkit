/*
  MeshTK.cc: to manage IO 

  created by bobye, Sep 11

  
 */

#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>





typedef TCLAP::ValueArg<std::string>                 Opt_String;
typedef TCLAP::ValueArg<double>                      Opt_scalar;
typedef TCLAP::ValueArg<int>                         Opt_Int;
typedef TCLAP::SwitchArg                             Opt_Bool;


int main(int argc, char** argv){

  try {
    TCLAP::CmdLine cmd("3D Mesh Toolkit, Supporting mesh (pre-)processing, visualization, ... \nContributor: Jianbo YE<yelpoo@gmail.com>", ' ', "0.1");

    Opt_String inputMeshName("i", "input_mesh", "File name of input mesh without file extension", true, "", "string", cmd); 

    Opt_String outputMeshName("o", "output_mesh", "File name of output mesh without file extension", false, "out", "string", cmd);

    Opt_String inputMeshType("", "input_mesh_type", "Input mesh file format, candidates are off(default), ...", false, "off", "string", cmd);

    Opt_String outputMeshType("", "output_mesh_type", "output mesh file format, candidates are off(default), ...", false, "off", "string", cmd);



    cmd.parse( argc, argv );















  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







