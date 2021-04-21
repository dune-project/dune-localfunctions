add_subdirectory(localfunctions)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
endif()
