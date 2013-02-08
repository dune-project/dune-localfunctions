set(TESTS test-biorthogonality)
add_executable(test-biorthogonality test-biorthogonality.cc)
target_link_libraries(test-biorthogonality dunegeometry dunecommon)
add_test(test-biorthogonality test-biorthogonality)

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
