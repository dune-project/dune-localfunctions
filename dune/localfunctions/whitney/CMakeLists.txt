set(whitneydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/whitney)
set(whitney_HEADERS
        edges0.5.hh)

set(SUBDIRS
        edges0.5)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${whitney_HEADERS} DESTINATION ${whitneydir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
