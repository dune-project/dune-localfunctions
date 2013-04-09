set(utilitydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/utility)

set(utility_HEADERS  basisevaluator.hh
                  basismatrix.hh
                  basisprint.hh
                  coeffmatrix.hh
                  defaultbasisfactory.hh
                  dglocalcoefficients.hh
                  field.hh
                  interpolationhelper.hh
                  l2interpolation.hh
                  lfematrix.hh
                  localfiniteelement.hh
                  monomialbasis.hh
                  multiindex.hh
                  polynomialbasis.hh
                  tensor.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${utility_HEADERS} DESTINATION ${utilitydir})
