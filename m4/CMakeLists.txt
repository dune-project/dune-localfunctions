set(ALLM4S  dune-localfunctions.m4)

set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/aclocal)
set(aclocal_DATA  ${ALLM4S})

set(EXTRA_DIST  ${ALLM4S})

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})
