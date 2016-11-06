// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGECOEFFICIENTS_HH
#define DUNE_LAGRANGECOEFFICIENTS_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  template< template <class,unsigned int> class LP, unsigned int dim, class F >
  struct LagrangeCoefficientsFactory;

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F>
  struct LagrangeCoefficientsFactoryTraits
  {
    static const unsigned int dimension = dim;
    const typedef LP<F,dim> Object;
    typedef unsigned int Key;
    typedef LagrangeCoefficientsFactory< LP,dim,F > Factory;
  };

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F>
  struct LagrangeCoefficientsFactory :
    public TopologyFactory< LagrangeCoefficientsFactoryTraits< LP,dim,F> >
  {
    typedef LagrangeCoefficientsFactoryTraits<LP,dim,F> Traits;
    static const unsigned int dimension = dim;
    typedef typename Traits::Object Object;
    typedef typename Traits::Key Key;

    template< class T >
    static Object *createObject ( const Key &order )
    {
      if (order == 0 || !Object::template supports<T>(order))
        return 0;
      typedef typename std::remove_const<Object>::type LagrangeCoefficients;
      LagrangeCoefficients *object = new LagrangeCoefficients(order);
      if ( !object->template build<T>() )
      {
        delete object;
        object = nullptr;
      }
      return object;
    }
  };

}

#endif // DUNE_LAGRANGECOEFFICIENTS_HH
