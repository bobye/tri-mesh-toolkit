#define _CRT_SECURE_NO_DEPRECATE

#if defined(_MSC_VER)
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR		char
# define _T(x)		x
# define _tprintf	printf
# define _tmain		main
#endif
#include <stdio.h>
#include <locale.h>

#include "simpleopt/SimpleOpt.h"

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FILE, OPT_VIEW, OPT_EXPL, OPT_LAPL, OPT_EMBD};

CSimpleOpt::SOption g_rgOptions[] = {
  { OPT_FILE,  _T("-f"),     SO_REQ_SEP }, // 
  { OPT_VIEW,  _T("-v"),     SO_NONE    },
  { OPT_EXPL,  _T("-e"),     SO_NONE    },
  { OPT_LAPL,  _T("-l"),     SO_NONE    },
  { OPT_EMBD,  _T("-b"),     SO_NONE    },  
  { OPT_HELP,  _T("-?"),     SO_NONE    }, // "-?"
  { OPT_HELP,  _T("--help"), SO_NONE    }, // "--help"
  SO_END_OF_OPTIONS                       // END
};

void ShowUsage() {
  _tprintf(_T("Usage: shape-caft [-f FILE] [-v] [-?] [--help] \n")
	   _T("\n")
	   _T("-f FILE_PREFIX      Load single shape file, binary/grayscale image\n")
	   _T("-v                  View labeled mesh\n")
	   _T("-e                  Load examples, assembly and export GraphCut matrix\n")
	   _T("-f                  Assembly and export linear LB matrix\n")
	   _T("-b                  Solve eigenvectors, and export embedding\n")
	   );
}

#include "meshtk/ManifoldTriMesh.hh"
#include "meshtk/MeshViewer.hh"

int main(int argc, char *argv[])
{
  CSimpleOpt args(argc, argv, g_rgOptions);
  meshtk::ManifoldTriMesh mesh;
  std::string nameprefix;
  bool view_label =false, load_examples = false, export_embed = false, assembly_LB = false;
  while (args.Next()) {
    if (args.LastError() == SO_SUCCESS) {

      switch (args.OptionId()) {
      case OPT_HELP:
	ShowUsage(); return 0;
      case OPT_FILE:	  
	nameprefix = std::string(args.OptionArg());
	//in_img = cv::imread(args.OptionArg(), 0);
	break;
      case OPT_VIEW:
	view_label = true; break;
      case OPT_EXPL:
	load_examples = true;	break;
      case OPT_LAPL:
	assembly_LB = true; break;
      case OPT_EMBD:
	export_embed = true; break;
      default:
	break;
      }
    }
    else {
      _tprintf(_T("Invalid argument: %s\n"), args.OptionText());
      return 1;
    }
  }

  std::string  reference_mesh_name = nameprefix, 
    deformation_mesh_name = nameprefix, 
    manifold_mesh_name = nameprefix,
    cluster_mesh_name = nameprefix,
    examples_mesh_name = nameprefix;
  if (nameprefix.compare("") != 0) {
    reference_mesh_name.append("_reference");
    deformation_mesh_name.append("_deformation");
    examples_mesh_name.append("_examples");

    mesh.read(reference_mesh_name, "off");
    mesh.init_index();// initialize mesh    
    mesh.update_base();
    mesh.update_compact_base();

    if (view_label) {
      cluster_mesh_name.append("_cluster");
      mesh.load_proxy_bone(cluster_mesh_name);

      meshtk::MeshViewer viewer(argc, argv);
      meshtk::MeshLabel painter(&mesh, mesh.vertex_label, mesh.label_color);
      //meshtk::MeshPainter painter(&mesh);

      viewer.add_painter(&painter);
      viewer.init();
      viewer.view();            
    } else if (load_examples) {
      mesh.load_examples(examples_mesh_name);
      mesh.PETSc_init(argc, argv);
      mesh.PETSc_assemble_graphcut();
      mesh.PETSc_export_LBmat(examples_mesh_name);
      mesh.PETSc_destroy();
    } else if (assembly_LB) {
      mesh.PETSc_init(argc, argv);
      mesh.PETSc_assemble_linearFEM_LBmat();
      mesh.PETSc_export_LBmat(examples_mesh_name);
      mesh.PETSc_destroy();
    }else if (export_embed) {
      mesh.PETSc_init(argc, argv);
      mesh.PETSc_load_LBmat(examples_mesh_name);
      mesh.PETSc_load_LBeigen(examples_mesh_name);
      mesh.PETSc_export_graphcut_vectors(examples_mesh_name);
      mesh.PETSc_destroy();      
    } else {
      mesh.load_sequence(deformation_mesh_name);
      mesh.compute_rotate_sequence();

      manifold_mesh_name.append("_manifold");
      mesh.print_rotate_sequence(manifold_mesh_name);
    }
  }


  
  return 0;
}






