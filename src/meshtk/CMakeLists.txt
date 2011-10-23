SET(MESHTK_LIB_SRC mesh_assist.cc TriMesh_Base.cc TriMesh_Dist.cc TriMesh_SIFT.cc DynamicTriMesh.cc MeshViewer.cc TriMesh_Curv.cc TriMesh_PETSc.cc TriMesh_UI.cc)

ADD_LIBRARY(meshtk SHARED ${MESHTK_LIB_SRC})
TARGET_LINK_LIBRARIES(meshtk ${LIBS})
