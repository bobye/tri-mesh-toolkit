#ifndef _MESH_IO_H
#define _MESH_IO_H
#include <CGAL/IO/Polyhedron_iostream.h>

struct Mesh_i{
  template <class Polyhedron>
  void operator()(Polyhedron &P, std::string file, std::string type){
    if (type.compare("off")==0){
      std::ifstream mesh_Fin;
      file.append("."); file.append("off");    
      mesh_Fin.open(file.c_str());
      if (!mesh_Fin.is_open())
	std::cerr << "Cannot open file " << file << std::endl;
      mesh_Fin >> P;// mesh read
      mesh_Fin.close();
    }
  };
};

struct Mesh_o{
  template <class Polyhedron>
  void operator()(Polyhedron &P, std::string file, std::string type){
    if (type.compare("off")==0){
      std::ofstream mesh_Fout;
      file.append("."); file.append(type);
      mesh_Fout.open(file.c_str());
      mesh_Fout << P;// mesh output
      std::cout << "Export mesh to: " << file <<  std::endl;
      mesh_Fout.close();
    }
  };
};

#endif
