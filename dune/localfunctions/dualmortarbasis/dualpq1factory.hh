// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_DUAL_P1_Q1_FACTORY_HH
#define DUNE_LOCALFUNCTIONS_DUAL_P1_Q1_FACTORY_HH

#include <map>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/dualmortarbasis.hh>

namespace Dune {

template<class D, class R, int dim, bool faceDual=false>
class DualPQ1LocalFiniteElementCache
{
protected:
  typedef Dune::DualP1LocalFiniteElement<D,R,dim,faceDual> DualP1;
  typedef Dune::DualQ1LocalFiniteElement<D,R,dim,faceDual> DualQ1;
  typedef typename DualP1::Traits::LocalBasisType::Traits T;
  typedef Dune::LocalFiniteElementVirtualInterface<T> FE;
  typedef std::map<Dune::GeometryType,FE*> FEMap;

public:
  /** \brief Type of the finite elements stored in this cache */
  typedef FE FiniteElementType;

  ~DualPQ1LocalFiniteElementCache()
  {
    typename FEMap::iterator it = cache_.begin();
    typename FEMap::iterator end = cache_.end();
    for(; it!=end; ++it)
      delete it->second;
  }

  //! create finite element for given GeometryType
  static FE* create(const Dune::GeometryType& gt)
  {
    if (gt.isSimplex())
      return new Dune::LocalFiniteElementVirtualImp<DualP1>(DualP1());
    if (gt.isCube())
      return new Dune::LocalFiniteElementVirtualImp<DualQ1>(DualQ1());
    return 0;
  }

  //! Get local finite element for given GeometryType
  const FiniteElementType& get(const Dune::GeometryType& gt) const
  {
    typename FEMap::const_iterator it = cache_.find(gt);
    if (it==cache_.end())
    {
      FiniteElementType* fe = create(gt);

      if (fe==0)
        DUNE_THROW(Dune::NotImplemented,"No Dual P/Q1 like local finite element available for geometry type " << gt);

      cache_[gt] = fe;
      return *fe;
    }
    return *(it->second);
  }

protected:
  mutable FEMap cache_;
};

}  // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_DUAL_P1_Q1_FACTORY_HH
