add_subdirectory(dualp1)
add_subdirectory(dualq1)
add_subdirectory(test EXCLUDE_FROM_ALL)

install(FILES
  dualp1.hh
  dualq1.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/dualmortarbasis)

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
