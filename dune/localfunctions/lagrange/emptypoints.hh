// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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

    Vector point_;
    LocalKey localKey_;
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
    EmptyPointSet ( const unsigned int order )
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

    unsigned int order () const
    {
      return order_;
    }

    unsigned int size () const
    {
      return points_.size();
    }

  protected:
    unsigned int order_;
    std::vector< LagrangePoint > points_;
  };

}

#endif // DUNE_LAGRANGE_EMPTYPOINTS_HH
