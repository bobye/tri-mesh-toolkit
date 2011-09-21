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

//#include "mesh_topo.hh"
#include "TriMesh.hh"
#include "MeshViewer.hh"

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
    /***************************************************************************/
    //contruct local chart for vertices and facets
    
    //mesh.update_vertex_localchart();
    //mesh.update_facet_localchart();
    /***************************************************************************/
    //Compute facet curvature tensor
    //mesh.update_facet_curvature();



    //mesh.write(output_mesh_name, output_mesh_type);

    MeshViewer viewer(argc, argv);
    MeshPainter painter(&mesh);
    viewer.add_painter(&painter);

    viewer.init();// call this func last before loop
    viewer.view();

  } catch (TCLAP::ArgException &e) //catch any expections 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }


  return 0;
}







