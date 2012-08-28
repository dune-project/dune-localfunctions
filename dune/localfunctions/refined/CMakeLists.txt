set(refineddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined)
set(refined_HEADERS  refinedp0.hh refinedp1.hh)

set(SUBDIRS  common refinedp0 refinedp1)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${refined_HEADERS} DESTINATION ${refineddir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
