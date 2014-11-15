set(monomdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/monom)

set(monom_HEADERS  monomlocalbasis.hh monomlocalcoefficients.hh monomlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${monom_HEADERS} DESTINATION ${monomdir})
