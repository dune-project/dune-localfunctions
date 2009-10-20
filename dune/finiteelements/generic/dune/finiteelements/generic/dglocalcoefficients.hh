// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGLOCALCOEFFICIENTS_HH
#define DUNE_DGLOCALCOEFFICIENTS_HH

#include <cassert>
#include <vector>

#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/finiteelements/common/localcoefficients.hh>

namespace Dune
{

  // DGLocalCoefficients
  // -------------------

  class DGLocalCoefficients
    : public LocalCoefficientsInterface< DGLocalCoefficients >
  {
    typedef DGLocalCoefficients This;
    typedef LocalCoefficientsInterface< This > Base;

  public:
    DGLocalCoefficients ( const unsigned int n )
      : localKey_( n )
    {
      for( unsigned i = 0; i < n; ++i )
        localKey_[ i ] = LocalKey( 0, 0, i );
    }

    const LocalKey &localKey ( const unsigned int i ) const
    {
      assert( i < size() );
      return localKey_[ i ];
    }

    unsigned int size () const
    {
      return localKey_.size();
    }

  private:
    std::vector< LocalKey > localKey_;
  };



  // DGLocalCoefficientsFactory
  // --------------------------
  template< class BasisCreator >
  struct DGLocalCoefficientsFactory;
  template< class BasisFactory >
  struct DGLocalCoefficientsFactoryTraits
  {
    static const unsigned int dimension = BasisFactory::dimension;
    typedef typename BasisFactory::Key Key;
    typedef DGLocalCoefficients LocalCoefficients;
    typedef const DGLocalCoefficients Object;
    typedef DGLocalCoefficientsFactory<BasisFactory> Factory;
  };

  template< class BasisFactory >
  struct DGLocalCoefficientsFactory :
    public TopologyFactory< DGLocalCoefficientsFactoryTraits<BasisFactory> >
  {
    typedef DGLocalCoefficientsFactoryTraits<BasisFactory> Traits;

    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      const typename BasisFactory::Basis *basis
        = BasisFactory::template create< Topology >( key );
      Object *coefficients = new Object( basis->size() );
      BasisFactory::release( basis );
      return coefficients;
    }
  };

}

#endif // #ifndef DUNE_DGLOCALCOEFFICIENTS_HH
