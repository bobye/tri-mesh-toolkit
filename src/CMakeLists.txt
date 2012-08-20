ADD_EXECUTABLE(MeshTK main.cc)
#ADD_EXECUTABLE(eigen_solver eigen_solver.cc)
TARGET_LINK_LIBRARIES(MeshTK meshtk ${LIBS})
#TARGET_LINK_LIBRARIES(eigen_solver ${SOLVER_LIBS})
