// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_PQKFACTORY_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_PQKFACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/p0.hh>

namespace Dune
{

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Empty default implementation
   */
  template<class D, class R, int d, int k>
  struct DimSpecificPQkLocalFiniteElementFactory
  {
    typedef typename P0LocalFiniteElement<D,R,d>::Traits::LocalBasisType::Traits T;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType&)
    {
      return nullptr;
    }
  };

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Specialization for dim=3
   */
  template<class D, class R, int k>
  struct DimSpecificPQkLocalFiniteElementFactory<D,R,3,k>
  {
    typedef typename P0LocalFiniteElement<D,R,3>::Traits::LocalBasisType::Traits T;
    using PrismP1 = LagrangePrismLocalFiniteElement<D,R,1>;
    using PrismP2 = LagrangePrismLocalFiniteElement<D,R,2>;
    using PyramidP1 = LagrangePyramidLocalFiniteElement<D,R,1>;
    using PyramidP2 = LagrangePyramidLocalFiniteElement<D,R,2>;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      if ((gt.isPrism())and (k==1))
        return new LocalFiniteElementVirtualImp<PrismP1>(PrismP1());
      if ((gt.isPrism())and (k==2))
        return new LocalFiniteElementVirtualImp<PrismP2>(PrismP2());
      if ((gt.isPyramid())and (k==1))
        return new LocalFiniteElementVirtualImp<PyramidP1>(PyramidP1());
      if ((gt.isPyramid())and (k==2))
        return new LocalFiniteElementVirtualImp<PyramidP2>(PyramidP2());
      return nullptr;
    }
  };


  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int k>
  struct PQkLocalFiniteElementFactory
  {
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    using P0 = P0LocalFiniteElement<D,R,dim>;
    using Pk = LagrangeSimplexLocalFiniteElement<D,R,dim,k>;
    using Qk = LagrangeCubeLocalFiniteElement<D,R,dim,k>;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      if (k==0)
        return new LocalFiniteElementVirtualImp<P0>(P0(gt));

      if (gt.isSimplex())
        return new LocalFiniteElementVirtualImp<Pk>(Pk());

      if (gt.isCube())
        return new LocalFiniteElementVirtualImp<Qk>(Qk());

      return DimSpecificPQkLocalFiniteElementFactory<D,R,dim,k>::create(gt);
    }
  };

}

#endif
