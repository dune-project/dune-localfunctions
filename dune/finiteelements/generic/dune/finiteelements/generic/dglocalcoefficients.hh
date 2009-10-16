// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGLOCALCOEFFICIENTS_HH
#define DUNE_DGLOCALCOEFFICIENTS_HH

#include <cassert>
#include <vector>

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



  // DGLocalCoefficientsCreator
  // --------------------------

  template< class BasisCreator >
  struct DGLocalCoefficientsCreator
  {
    static const unsigned int dimension = BasisCreator::dimension;

    typedef typename BasisCreator::Key Key;

    typedef DGLocalCoefficients LocalCoefficients;

    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &key )
    {
      const typename BasisCreator::Basis &basis
        = BasisCreator::template basis< Topology >( key );
      const LocalCoefficients *coefficients = new LocalCoefficients( basis.size() );
      BasisCreator::release( basis );
      return *coefficients;
    }

    static void release ( const LocalCoefficients &localCoefficients )
    {
      delete &localCoefficients;
    }

    template< class Topology >
    static bool supports ( const Key &key )
    {
      return true;
    }
  };

}

#endif // #ifndef DUNE_DGLOCALCOEFFICIENTS_HH
