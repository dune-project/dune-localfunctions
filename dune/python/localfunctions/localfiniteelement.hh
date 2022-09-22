// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_LOCALFUNCTIONS_LOCALFINITEELEMENT_HH
#define DUNE_PYTHON_LOCALFUNCTIONS_LOCALFINITEELEMENT_HH

#include <dune/python/pybind11/pybind11.h>

#include <dune/common/visibility.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/virtualinterface.hh>

namespace Dune {
namespace Python {

namespace detail {

template<typename LocalBasis>
DUNE_EXPORT auto registerLocalBasis(pybind11::handle scope)
{
  static auto cls = pybind11::class_<LocalBasis>(scope, "LocalBasis");

  cls.def("__len__", [](const LocalBasis& basis) { return basis.size(); });
  cls.def_property_readonly("order", [](const LocalBasis& basis) { return basis.order(); });
  cls.def("evaluateFunction",
          [](const LocalBasis& basis, const typename LocalBasis::Traits::DomainType& in) {
            std::vector<typename LocalBasis::Traits::RangeType> out;
            basis.evaluateFunction(in, out);
            return out;
          });
  cls.def("evaluateJacobian",
          [](const LocalBasis& basis, const typename LocalBasis::Traits::DomainType& in) {
            std::vector<typename LocalBasis::Traits::JacobianType> out;
            basis.evaluateJacobian(in, out);
            return out;
          });
  return cls;
}

DUNE_EXPORT auto registerLocalKey(pybind11::handle scope)
{
  static auto cls = pybind11::class_<LocalKey>(scope, "LocalKey");

  cls.def_property_readonly("subEntity", &LocalKey::subEntity);
  cls.def_property_readonly("codim", &LocalKey::codim);
  cls.def_property("index",
                   [](const LocalKey& key) { return key.index(); },
                   [](LocalKey& key, unsigned int index) { key.index(index); });
  cls.def("__lt__", &LocalKey::operator<);

  return cls;
}

} /* namespace detail */

template<typename LocalFiniteElement>
DUNE_EXPORT auto registerLocalFiniteElement(pybind11::handle scope, const char* name = "LocalFiniteElement")
{
  static auto cls = pybind11::class_<LocalFiniteElement>(scope, name);

  detail::registerLocalBasis<typename LocalFiniteElement::Traits::LocalBasisType>(cls);

  cls.def_property_readonly("localBasis", &LocalFiniteElement::localBasis, pybind11::return_value_policy::reference_internal);
  // cls.def_property_readonly("localCoefficients", &LocalFiniteElement::localCoefficients);
  // cls.def_property_readonly("localInterpolation", &LocalFiniteElement::localInterpolation);
  cls.def("__len__", &LocalFiniteElement::size);
  cls.def_property_readonly("type", &LocalFiniteElement::type);

  return cls;
}


} /* namespace Python */
} /* namespace Dune */

#endif
