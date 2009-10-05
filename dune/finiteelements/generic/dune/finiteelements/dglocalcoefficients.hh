// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGLOCALCOEFFICIENTS_HH
#define DUNE_DGLOCALCOEFFICIENTS_HH

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

  template< class NumBasisFunctions >
  struct DGLocalCoefficientsCreator
  {
    static const unsigned int dimension = NumBasisFunctions::dimension;

    typedef typename NumBasisFunctions::Key Key;

    typedef DGLocalCoefficients LocalCoefficients;

    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &key )
    {
      const unsigned int numCoefficients
        = NumBasisFunctions::template numBasisFunctions< Topology >( key );
      return *(new LocalCoefficients( numCoefficients ));
    }

    static void release ( const LocalCoefficients &localCoefficients )
    {
      delete &localCoefficients;
    }
  };

}

#endif // #ifndef DUNE_DGLOCALCOEFFICIENTS_HH
