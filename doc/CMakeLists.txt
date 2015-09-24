add_subdirectory(doxygen)

dune_add_latex_document(dune-localfunctions-manual.tex
  FATHER_TARGET doc
  DEFAULT_PDF)

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/dune-localfunctions-manual.pdf
  ${CMAKE_INSTALL_DOCDIR})
