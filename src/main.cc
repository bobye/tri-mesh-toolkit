/*
  MeshTK.cc: to manage IO 

  created by bobye, Sep 11

  
 */

#include <string>
#include <algorithm>

#include <tclap/CmdLine.h>
typedef TCLAP::ValueArg<std::string>                 Opt_String;
typedef TCLAP::ValueArg<double>                      Opt_scalar;
typedef TCLAP::ValueArg<int>                         Opt_Int;
typedef TCLAP::SwitchArg                             Opt_Bool;

//#include "mesh_topo.h"
#include "TriMesh.h"


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

    
    TriMesh mesh;// index mapping from array to halfedge data structure

    mesh.read(input_mesh_name, input_mesh_type);
    mesh.init_index();
    // Update mesh infomation
    mesh.update_base();

    //mesh processing
    Vec_Fun vertex_LC[2] = {Vec_Fun(mesh.vertex_num), Vec_Fun(mesh.vertex_num)};
    Vec_Fun facet_LC[2]={Vec_Fun(mesh.facet_num), Vec_Fun(mesh.facet_num)};
    

    mesh.init_vertex_localchart(vertex_LC[0], vertex_LC[1]); 
    mesh.init_facet_localchart(facet_LC[0], facet_LC[1]);


    
    mesh.write(output_mesh_name, output_mesh_type);
  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







