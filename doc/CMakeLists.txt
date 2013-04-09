# $Id: Makefile.am 836 2007-12-05 17:34:43Z robertk $

set(SUBDIRS  doxygen)

set(TEXSOURCES  dune-localfunctions-manual.tex)

if( BUILD_DOCS)
  set(DOCFILES  dune-localfunctions-manual.pdf)
  set(EXTRA_DIST  ${DOCFILES})
  set(EXTRAINSTALL  ${DOCFILES})
endif( BUILD_DOCS )

# setting like in dune-web
set(BASEDIR ../..)
set(CURDIR doc/localfunctions)
set(docdir ${CMAKE_INSTALL_DOCDIR}/dune-localfunctions)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
# include not needed for CMake
# include $(top_srcdir)/am/webstuff
# include not needed for CMake
# include $(top_srcdir)/am/latex

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

# Install documentation
# We assume that ${DOCFILES} and ${DOCFILES_EXTRA} are targets
install(TARGETS ${DOCFILES} ${DOCFILES_EXTRA} DESTINATION ${docdir})
