// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH

#include <sstream>

#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/generic/basisprint.hh>
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

    typedef Dune::MonomialBasisProvider< dimension, StorageField > MonomialBasisProviderType;
    typedef typename MonomialBasisProviderType::Basis MonomialBasisType;

    typedef unsigned int Key;

    // typedef amp::ampf< Precision< CF >::value > ComputeField;
    typedef CF ComputeField;
    typedef SparseCoeffMatrix< StorageField, 1 > CoefficientMatrix;
    typedef StandardEvaluator< MonomialBasisType > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;

    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    template< class Topology >
    static const Basis &basis ( const unsigned int order )
    {
      const MonomialBasisType &monomialBasis = MonomialBasisProviderType::template basis< SimplexTopology >( order );

      static CoefficientMatrix _coeffs;
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
