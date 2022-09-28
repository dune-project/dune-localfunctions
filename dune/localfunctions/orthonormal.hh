// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_ORTHONORMALFINITEELEMENT_HH
#define DUNE_ORTHONORMALFINITEELEMENT_HH

#include <dune/localfunctions/utility/localfiniteelement.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>
#include <dune/localfunctions/utility/l2interpolation.hh>
#include <dune/localfunctions/orthonormal/orthonormalbasis.hh>

namespace Dune
{
  /**
   * \brief A class providing orthonormal basis functions
   *
   * These basis functions are constructed by l2 orthonormalizing
   * of monomials over a reference element. Thus the span is
   * always Pk. The coefficients and the interpolation are given
   * by the Dune::DGLocalCoefficientsFactory and the
   * Dune::LocalL2InterpolationFactory.
   *
   * \ingroup Orthonormal
   *
   * \tparam dimDomain dimension of reference elements
   * \tparam D domain for basis functions
   * \tparam R range for basis functions
   * \tparam SF storage field for basis matrix
   * \tparam CF compute field for basis matrix
   **/
  template< unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class OrthonormalLocalFiniteElement
    : public GenericLocalFiniteElement< OrthonormalBasisFactory< dimDomain, SF, CF >,
          DGLocalCoefficientsFactory< OrthonormalBasisFactory< dimDomain, SF, CF > >,
          LocalL2InterpolationFactory< OrthonormalBasisFactory< dimDomain, SF, CF >,true > >
  {
    typedef GenericLocalFiniteElement< OrthonormalBasisFactory< dimDomain, SF, CF >,
        DGLocalCoefficientsFactory< OrthonormalBasisFactory< dimDomain, SF, CF > >,
        LocalL2InterpolationFactory< OrthonormalBasisFactory< dimDomain, SF, CF >,true > > Base;
  public:
    using typename Base::Traits;

    /** \todo Please doc me !
     */
    OrthonormalLocalFiniteElement ( const GeometryType &gt, unsigned int order )
      : Base(gt, order)
    {}
  };

}

#endif
