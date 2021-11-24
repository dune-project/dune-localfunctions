add_python_targets(localfunctions
  __init__
)

dune_add_pybind11_module(NAME _localfunctions)
set_property(TARGET _localfunctions PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
install(TARGETS _localfunctions LIBRARY DESTINATION python/dune/localfunctions)
