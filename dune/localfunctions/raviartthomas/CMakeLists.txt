# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(raviartthomas02d)
add_subdirectory(raviartthomas03d)
add_subdirectory(raviartthomas0cube2d)
add_subdirectory(raviartthomas0cube3d)
add_subdirectory(raviartthomas12d)
add_subdirectory(raviartthomas1cube2d)
add_subdirectory(raviartthomas1cube3d)
add_subdirectory(raviartthomas2cube2d)
add_subdirectory(raviartthomas3cube2d)
add_subdirectory(raviartthomas4cube2d)
add_subdirectory(raviartthomassimplex)

install(FILES
  raviartthomas02d.hh
  raviartthomas03d.hh
  raviartthomas0cube2d.hh
  raviartthomas0cube3d.hh
  raviartthomas12d.hh
  raviartthomas1cube2d.hh
  raviartthomas1cube3d.hh
  raviartthomas2cube2d.hh
  raviartthomas3cube2d.hh
  raviartthomas4cube2d.hh
  raviartthomascube.hh
  raviartthomaslfecache.hh
  raviartthomassimplex.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas)
