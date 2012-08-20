SET(MESHTK_LIB_SRC mesh_assist.cc   TriMesh_Curv.cc   TriMesh_SIFT.cc
  DynamicTriMesh.cc  MeshViewer.cc    TriMesh_Dist.cc   TriMesh_UI.cc
  lapack_wrapper.c   TriMesh_Base.cc  TriMesh_PETSc.cc )

ADD_LIBRARY(meshtk STATIC ${MESHTK_LIB_SRC})
TARGET_LINK_LIBRARIES(meshtk ${LIBS})
