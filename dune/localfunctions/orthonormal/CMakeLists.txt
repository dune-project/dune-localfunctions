set(orthonormaldir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/orthonormal)
set(orthonormal_HEADERS  orthonormalbasis.hh orthonormalcompute.hh)

install(FILES ${orthonormal_HEADERS} DESTINATION ${orthonormaldir})
