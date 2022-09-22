// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/localfunctions/localfiniteelement.hh>

PYBIND11_MODULE(_localfunctions, module)
{
  Dune::Python::detail::registerLocalKey(module);
}
