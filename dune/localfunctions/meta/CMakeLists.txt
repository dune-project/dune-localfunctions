set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/meta)
set(my_HEADERS
        power.hh)

set(SUBDIRS
        power)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${my_HEADERS} DESTINATION ${mydir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
