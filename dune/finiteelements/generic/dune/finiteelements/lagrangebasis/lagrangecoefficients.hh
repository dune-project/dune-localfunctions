// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGECOEFFICIENTS_HH
#define DUNE_LAGRANGECOEFFICIENTS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/finiteelements/common/field.hh>
#include <dune/common/forloop.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

#include <dune/finiteelements/common/localcoefficients.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< template <class,unsigned int> class LC, unsigned int dim, class F >
  class LagrangeCoefficientsFactory;


  // LagrangeCoefficient
  // -------------

  template< class F, unsigned int dim >
  class LagrangeCoefficient
  {
    typedef LagrangeCoefficient< F, dim > This;

    template< class, class >
    friend class LagrangeCoefficientsImpl;

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

    // private:
    Vector point_;
    LocalKey localKey_;
  };

  // LagnragePoints
  // --------------

  template< class F, unsigned int dim >
  class LagrangeCoefficientsBase
    : public LocalCoefficientsInterface< LagrangeCoefficientsBase< F, dim > >
  {
    typedef LagrangeCoefficientsBase< F, dim > This;

  public:
    typedef F Field;

    static const unsigned int dimension = dim;

    typedef Dune::LagrangeCoefficient< Field, dimension > LagrangeCoefficient;

    typedef typename std::vector< LagrangeCoefficient >::const_iterator iterator;

    // private:
    LagrangeCoefficientsBase ( const unsigned int order )
      : order_( order )
    {}

    void resize( const unsigned int size )
    {
      points_.resize(size);
    }

  public:
    const LagrangeCoefficient &operator[] ( const unsigned int i ) const
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

    // private:
    unsigned int order_;
    std::vector< LagrangeCoefficient > points_;
  };

  // LagrangeCoefficientsFactory
  // ---------------------
  template< template <class,unsigned int> class LC,
      unsigned int dim, class F>
  struct LagrangeCoefficientsFactoryTraits
  {
    static const unsigned int dimension = dim;
    const typedef LC<F,dim> Object;
    typedef unsigned int Key;
    typedef LagrangeCoefficientsFactory< LC,dim,F > Factory;
  };

  template< template <class,unsigned int> class LC,
      unsigned int dim, class F>
  struct LagrangeCoefficientsFactory :
    public TopologyFactory< LagrangeCoefficientsFactoryTraits< LC,dim,F> >
  {
    typedef LagrangeCoefficientsFactoryTraits<LC,dim,F> Traits;
    static const unsigned int dimension = dim;
    typedef typename Traits::Object Object;
    typedef typename Traits::Key Key;

    template< class T >
    static Object *createObject ( const Key &order )
    {
      if (!Object::template supports<T>(order))
        return 0;
      typedef typename remove_const<Object>::type LagrangeCoefficients;
      LagrangeCoefficients *object = new LagrangeCoefficients(order);
      if ( !object->template build<T>() )
      {
        delete object;
        object = 0;
      }
      return object;
    }
  };

}

#endif
