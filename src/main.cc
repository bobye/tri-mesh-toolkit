/*
  MeshTK.cc: to manage IO 

  created by bobye, Sep 11

  
 */

#include <string>
#include <iostream>
#include <fstream>

#include <algorithm>

#include <tclap/CmdLine.h>
typedef TCLAP::ValueArg<std::string>                 Opt_String;
typedef TCLAP::ValueArg<double>                      Opt_scalar;
typedef TCLAP::ValueArg<int>                         Opt_Int;
typedef TCLAP::SwitchArg                             Opt_Bool;

#include <CGAL/IO/Polyhedron_iostream.h>
#include "mesh_topo.h"

Polyhedron_IS PI;// index mapping from array to halfedge data structure



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

    std::ifstream mesh_Fin;
    
    std::string file_name_in(input_mesh_name);
    file_name_in.append("."); file_name_in.append(input_mesh_type);
    
    mesh_Fin.open(file_name_in.c_str());
    if (!mesh_Fin.is_open())
      std::cerr << "Cannot open file " << file_name_in.c_str() << std::endl;
    mesh_Fin >> PI.P;// mesh read
    mesh_Fin.close();
    Polyhedron_Init()(PI);// init all data of PI, with P


    //mesh processing



    


    std::ofstream mesh_Fout;
    std::string file_name_out(output_mesh_name);
    file_name_out.append("."); file_name_out.append(output_mesh_type);
    mesh_Fout.open(file_name_out.c_str());
    mesh_Fout << PI.P;// mesh output
    std::cout << "Put output mesh: " << file_name_out.c_str() <<  std::endl;
    mesh_Fout.close();
  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







