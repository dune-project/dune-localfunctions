set(prismp2dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/prismp2)

set(prismp2_HEADERS  prismp2localbasis.hh prismp2localcoefficients.hh prismp2localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${prismp2_HEADERS} DESTINATION ${prismp2dir})
