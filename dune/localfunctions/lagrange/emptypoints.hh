// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LAGRANGE_EMPTYPOINTS_HH
#define DUNE_LAGRANGE_EMPTYPOINTS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  template< class F, unsigned int dim >
  class LagrangePoint
  {
    typedef LagrangePoint< F, dim > This;

    template< class, class >
    friend class LagrangePointSetImpl;

  public:
    static const int dimension = dim;

    typedef F Field;

    typedef FieldVector< Field, dimension > Vector;

    const Vector &point () const
    {
      return point_;
    }

    const LocalKey &localKey () const
    {
      return localKey_;
    }

    const Field weight () const
    {
      return weight_;
    }

    Vector point_ = {};
    LocalKey localKey_ = {};
    Field weight_ = {};
  };

  // EmptyPointSet
  // --------------

  template< class F, unsigned int dim >
  class EmptyPointSet
  {
    typedef EmptyPointSet< F, dim > This;

  public:
    typedef F Field;

    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoint< Field, dimension > LagrangePoint;

    typedef typename std::vector< LagrangePoint >::const_iterator iterator;

  protected:
    EmptyPointSet ( const std::size_t order )
      : order_( order )
    {}

  public:
    const LagrangePoint &operator[] ( const unsigned int i ) const
    {
      assert( i < size() );
      return points_[ i ];
    }

    iterator begin () const
    {
      return points_.begin();
    }

    iterator end () const
    {
      return points_.end();
    }

    const LocalKey &localKey ( const unsigned int i ) const
    {
      return (*this)[ i ].localKey();
    }

    std::size_t order () const
    {
      return order_;
    }

    std::size_t size () const
    {
      return points_.size();
    }

  protected:
    std::size_t order_;
    std::vector< LagrangePoint > points_;
  };

}

#endif // DUNE_LAGRANGE_EMPTYPOINTS_HH
