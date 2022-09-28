# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(doxygen)

dune_add_latex_document(
  SOURCE dune-localfunctions-manual.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR}/comm)
