add_subdirectory(edges0.5)

install(FILES
  edges0.5.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/whitney)
