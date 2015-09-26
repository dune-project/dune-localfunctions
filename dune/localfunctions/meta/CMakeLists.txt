add_subdirectory(power)

install(FILES
  power.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/meta)
