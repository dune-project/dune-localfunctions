// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH

#include <sstream>

#include <dune/finiteelements/generic/topologyfactory.hh>

#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/generic/basisprint.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalcompute.hh>

namespace Dune
{

  // OrthonormalBasisFactory
  // -----------------------
  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisFactory;
  template< int dim, class SF, class CF >
  struct OrthonormalBasisFactoryTraits
  {
    typedef Dune::MonomialBasisProvider< dim, SF > MonomialBasisProviderType;
    typedef typename MonomialBasisProviderType::Object MonomialBasisType;
    typedef SparseCoeffMatrix< SF, 1 > CoefficientMatrix;
    typedef StandardEvaluator< MonomialBasisType > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;

    static const unsigned int dimension = dim;
    typedef unsigned int Key;
    typedef const Basis Object;
    typedef OrthonormalBasisFactory<dim,SF,CF> Factory;
  };

  template< int dim, class SF, class CF >
  struct OrthonormalBasisFactory :
    public TopologyFactory< OrthonormalBasisFactoryTraits<dim,SF,CF> >
  {
    static const unsigned int dimension = dim;
    typedef SF StorageField;
    typedef CF ComputeField;
    typedef OrthonormalBasisFactoryTraits<dim,SF,CF> Traits;

    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;

    template <unsigned int dd, class FF>
    struct EvaluationBasisFactory
    {
      typedef MonomialBasisProvider<dd,FF> Type;
    };

    typedef typename EvaluationBasisFactory< dimension, StorageField >::Type MonomialBasisProviderType;
    typedef typename MonomialBasisProviderType::Object MonomialBasisType;

    typedef SparseCoeffMatrix< StorageField, 1 > CoefficientMatrix;
    typedef StandardEvaluator< MonomialBasisType > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;

    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    template< class Topology >
    static Object *createObject ( const unsigned int order )
    {
      const typename Traits::MonomialBasisType &monomialBasis = *Traits::MonomialBasisProviderType::template create< SimplexTopology >( order );

      static typename Traits::CoefficientMatrix _coeffs;
      if( _coeffs.size() <= monomialBasis.size() )
      {
        ONBCompute::ONBMatrix< Topology, ComputeField > matrix( order );
        _coeffs.fill( matrix );
#if GLFEM_BASIS_PRINT
        {
          typedef MultiIndex< dimension > MIField;
          typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
          typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,dimension> > BasisMI;
          const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<Topology>(order+1);
          BasisMI basisMI(_mBasisMI);
          basisMI.fill(matrix);
          std::stringstream name;
          name << "orthonormal_" << Topology::name() << "_p" << order << ".basis";
          std::ofstream out(name.str().c_str());
          basisPrint<0>(out,basisMI);
        }
#endif
      }

      return new Basis( monomialBasis, _coeffs, monomialBasis.size() );
    }
  };

}

#endif // #ifndef DUNE_ORTHONORMALBASIS_HH
