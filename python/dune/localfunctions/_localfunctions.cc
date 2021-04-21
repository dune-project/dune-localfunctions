#include <dune/python/pybind11/pybind11.h>
#include <dune/python/localfunctions/localfiniteelement.hh>

PYBIND11_MODULE(_localfunctions, module)
{
  Dune::Python::detail::registerLocalKey(module);
}
