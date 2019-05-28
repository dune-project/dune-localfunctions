add_subdirectory(doxygen)

dune_add_latex_document(
  SOURCE dune-localfunctions-manual.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR}/comm)

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/dune-localfunctions-manual.pdf
  ${CMAKE_INSTALL_DOCDIR} dune-localfunctions-manual_safepdf)
