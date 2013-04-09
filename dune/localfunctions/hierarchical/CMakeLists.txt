set(hierarchicaldir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/hierarchical)
set(hierarchical_HEADERS  hierarchicalp2.hh hierarchicalp2withelementbubble.hh hierarchicalprismp2.hh)

set(SUBDIRS  hierarchicalp2 hierarchicalprismp2 hierarchicalp2withelementbubble)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${hierarchical_HEADERS} DESTINATION ${hierarchicaldir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
