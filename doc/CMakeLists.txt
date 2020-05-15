add_subdirectory(doxygen)

dune_add_latex_document(
  SOURCE dune-localfunctions-manual.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR}/comm)
