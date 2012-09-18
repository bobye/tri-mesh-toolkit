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
enum { OPT_HELP, OPT_FILE, OPT_VIEW};

CSimpleOpt::SOption g_rgOptions[] = {
  { OPT_FILE,  _T("-f"),     SO_REQ_SEP }, // 
  { OPT_VIEW,  _T("-v"),     SO_NONE    },
  { OPT_HELP,  _T("-?"),     SO_NONE    }, // "-?"
  { OPT_HELP,  _T("--help"), SO_NONE    }, // "--help"
  SO_END_OF_OPTIONS                       // END
};

void ShowUsage() {
  _tprintf(_T("Usage: shape-caft [-f FILE] [-v] [-?] [--help] \n")
	   _T("\n")
	   _T("-f FILE_PREFIX      Load single shape file, binary/grayscale image\n")
	   _T("-v                  View the transformed shape")
	   );
}

#include "meshtk/ManifoldTriMesh.hh"

int main(int argc, char *argv[])
{
  CSimpleOpt args(argc, argv, g_rgOptions);
  meshtk::ManifoldTriMesh mesh;
  std::string nameprefix, reference_mesh_name, deformation_mesh_name;
    while (args.Next()) {
        if (args.LastError() == SO_SUCCESS) {

	    switch (args.OptionId()) {
	    case OPT_HELP:
	      ShowUsage(); return 0;
	    case OPT_FILE:	  
	      nameprefix = std::string(args.OptionArg());
	      reference_mesh_name = nameprefix; 
	      deformation_mesh_name = nameprefix; 
	      reference_mesh_name.append("_reference");
	      deformation_mesh_name.append("_deformation");

	      mesh.read(reference_mesh_name, "off");
	      mesh.init_index();// initialize mesh    
	      mesh.update_base();
	      mesh.load_sequence(deformation_mesh_name);

	      //in_img = cv::imread(args.OptionArg(), 0);
	      break;
	    case OPT_VIEW:
	      break;
	    default:
	      break;
	    }
        }
        else {
	  _tprintf(_T("Invalid argument: %s\n"), args.OptionText());
	  return 1;
        }
    }


  
  return 0;
}
