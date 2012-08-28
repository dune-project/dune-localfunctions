set(brezzidouglasmarini1q2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1q2d)

set(brezzidouglasmarini1q2d_HEADERS        brezzidouglasmarini1q2dlocalbasis.hh
                                        brezzidouglasmarini1q2dlocalcoefficients.hh
                                        brezzidouglasmarini1q2dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${brezzidouglasmarini1q2d_HEADERS} DESTINATION ${brezzidouglasmarini1q2ddir})
