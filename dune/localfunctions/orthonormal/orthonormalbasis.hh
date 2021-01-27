// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH

#include <sstream>

#include <dune/localfunctions/utility/polynomialbasis.hh>
#include <dune/localfunctions/orthonormal/orthonormalcompute.hh>

namespace Dune
{

  // OrthonormalBasisFactory
  // -----------------------
  template< int dim, class D, class R, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisFactory
  {
    static const unsigned int dimension = dim;
    typedef D  Domain;
    typedef R  Range;
    typedef SF StorageField;
    typedef CF ComputeField;

    template <unsigned int dd, class FF>
    struct EvaluationBasisFactory
    {
      typedef MonomialBasisProvider<dd,FF> Type;
    };

    typedef typename EvaluationBasisFactory< dimension, StorageField >::Type MonomialBasisProviderType;
    typedef typename MonomialBasisProviderType::Object MonomialBasisType;

    typedef SparseCoeffMatrix< StorageField, 1 > CoefficientMatrix;
    typedef StandardEvaluator< MonomialBasisType > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix, Domain, Range > Basis;

    typedef unsigned int Key;
    typedef const Basis Object;

    typedef typename Impl::SimplexTopology< dim >::type SimplexTopology;

    template< class Topology >
    static Object *create ( const unsigned int order )
    {
      const MonomialBasisType &monomialBasis = *MonomialBasisProviderType::template create< SimplexTopology >( order );

      static CoefficientMatrix _coeffs;
      if( _coeffs.size() <= monomialBasis.size() )
      {
        ONBCompute::ONBMatrix< Topology, ComputeField > matrix( order );
        _coeffs.fill( matrix );
      }

      return new Basis( monomialBasis, _coeffs, monomialBasis.size() );
    }
    static void release( Object *object ) { delete object; }
  };

}

#endif // #ifndef DUNE_ORTHONORMALBASIS_HH
