set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/common)
set(common_HEADERS
        interface.hh
        interfaceswitch.hh
        localbasis.hh
        localkey.hh
        localfiniteelementtraits.hh
        localtoglobaladaptors.hh
        virtualinterface.hh
        virtualwrappers.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})
