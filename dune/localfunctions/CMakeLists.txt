# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(brezzidouglasfortinmarini)
add_subdirectory(brezzidouglasmarini)
add_subdirectory(common)
add_subdirectory(dualmortarbasis)
add_subdirectory(hierarchical)
add_subdirectory(lagrange)
add_subdirectory(meta)
add_subdirectory(mimetic)
add_subdirectory(monomial)
add_subdirectory(nedelec)
add_subdirectory(orthonormal)
add_subdirectory(rannacherturek)
add_subdirectory(raviartthomas)
add_subdirectory(refined)
add_subdirectory(test)
add_subdirectory(utility)
add_subdirectory(whitney)

install(FILES
  brezzidouglasmarini.hh
  crouzeixraviart.hh
  dualmortarbasis.hh
  hierarchical.hh
  lagrange.hh
  mimetic.hh
  monomial.hh
  nedelec.hh
  orthonormal.hh
  rannacherturek.hh
  raviartthomas.hh
  refined.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions)

# Install some test headers, because they might be useful for tests in
# other modules.
# We do this here as test will not be considered for make install.
install(FILES
  test/geometries.hh
  test/test-fe.hh
  test/test-localfe.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/test)
