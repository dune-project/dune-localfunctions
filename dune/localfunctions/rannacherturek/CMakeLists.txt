# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(rannacherturek2d)
add_subdirectory(rannacherturek3d)

install(FILES ${rannacherturek_HEADERS}
  rannacherturek.hh
  rannachertureklocalbasis.hh
  rannachertureklocalcoefficients.hh
  rannachertureklocalinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/rannacherturek)
