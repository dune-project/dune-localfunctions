set(brezzidouglasmarini12ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini12d)

set(brezzidouglasmarini12d_HEADERS  brezzidouglasmarini12dlocalbasis.hh
                                 brezzidouglasmarini12dlocalcoefficients.hh
                                 brezzidouglasmarini12dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${brezzidouglasmarini12d_HEADERS} DESTINATION ${brezzidouglasmarini12ddir})
