set(TESTS test-biorthogonality)
add_executable(test-biorthogonality test-biorthogonality.cc)
target_link_libraries(test-biorthogonality dunegeometry dunecommon)
add_test(test-biorthogonality test-biorthogonality)
