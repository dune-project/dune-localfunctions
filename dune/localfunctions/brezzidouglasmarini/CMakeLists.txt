# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(brezzidouglasmarini1cube2d)
add_subdirectory(brezzidouglasmarini1cube3d)
add_subdirectory(brezzidouglasmarini2cube2d)
add_subdirectory(brezzidouglasmarini1simplex2d)
add_subdirectory(brezzidouglasmarini2simplex2d)

install(FILES
  brezzidouglasmarinicube.hh
  brezzidouglasmarinisimplex.hh
  brezzidouglasmarini1cube2d.hh
  brezzidouglasmarini1cube3d.hh
  brezzidouglasmarini2cube2d.hh
  brezzidouglasmarini1simplex2d.hh
  brezzidouglasmarini2simplex2d.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini)
