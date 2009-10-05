// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH

#include <sstream>

#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalcompute.hh>

namespace Dune
{

  // OrthonormalBasisCreator
  // -----------------------

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisCreator
  {
    static const unsigned int dimension = dim;
    typedef SF StorageField;

    typedef Dune::MonomialBasisProvider< dimension, StorageField > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Basis MonomialBasis;

    typedef unsigned int Key;

    typedef AlgLib::MultiPrecision< Precision< CF >::value > ComputeField;
    typedef SparseCoeffMatrix< StorageField, 1 > CoefficientMatrix;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;

    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    template< class Topology >
    static const Basis &basis ( const unsigned int order )
    {
      const MonomialBasis &monomialBasis = MonomialBasisProvider::template basis< SimplexTopology >( order );

      static CoefficientMatrix _coeffs;
      if( _coeffs.size() <= monomialBasis.size() )
      {
        ONBCompute::ONBMatrix< Topology, ComputeField > matrix( order );
        _coeffs.fill( matrix );
      }

      return *(new Basis( monomialBasis, _coeffs, monomialBasis.size() ));
    }

    static void release ( const Basis &basis )
    {
      delete &basis;
    }
  };



  // OrthonormalBasisProvider
  // ------------------------

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisProvider
    : public BasisProvider< OrthonormalBasisCreator< dim, SF, CF > >
  {};

}

#endif // #ifndef DUNE_ORTHONORMALBASIS_HH
