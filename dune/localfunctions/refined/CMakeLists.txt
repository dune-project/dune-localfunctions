add_subdirectory(common)
add_subdirectory(refinedp0)
add_subdirectory(refinedp1)

install(FILES
  refinedp0.hh
  refinedp1.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined)
