// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/finiteelements/common/matrix.hh>
#include <dune/finiteelements/generic/basismatrix.hh>

#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/generic/multiindex.hh>
#include <dune/finiteelements/generic/monomialbasis.hh>
#include <dune/finiteelements/generic/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>

namespace Dune
{

  // LagrangeBasisFactory
  // --------------------

  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisFactory;
  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisTraits
  {
    typedef Dune::MonomialBasisProvider< dim, SF > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Object MonomialBasis;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< SF, 1 > > Basis;

    typedef Dune::MonomialBasisProvider< dim, CF > PreBasisFactory;
    typedef typename PreBasisFactory::Object PreBasis;
    typedef LagrangeInterpolationFactory< LP, dim, CF > InterpolationFactory;
    typedef typename InterpolationFactory::Object Interpolation;

    static const unsigned int dimension = dim;
    typedef const Basis Object;
    typedef unsigned int Key;
    typedef LagrangeBasisFactory<LP,dim,SF,CF> Factory;
  };

  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisFactory
    : public TopologyFactory< LagrangeBasisTraits< LP,dim,SF,CF > >
  {
    typedef LagrangeBasisTraits<LP,dim,SF,CF> Traits;
    static const unsigned int dimension = dim;
    typedef SF StorageField;
    typedef CF ComputeField;
    typedef typename Traits::Basis Basis;

    typedef typename Traits::Object Object;
    typedef typename Traits::Key Key;

    template< class Topology >
    static Object *createObject ( const Key &order )
    {
      const typename Traits::MonomialBasis *monomialBasis = Traits::MonomialBasisProvider::template create< Topology >( order );
      Basis *basis = new Basis( *monomialBasis );

      const typename Traits::PreBasis *preBasis = Traits::PreBasisFactory::template create<Topology>( order );
      const typename Traits::Interpolation *interpol = Traits::InterpolationFactory::template create<Topology>( order );
      BasisMatrix< typename Traits::PreBasis,
          typename Traits::Interpolation,
          ComputeField > matrix( *preBasis, *interpol );
      basis->fill( matrix );
      Traits::PreBasisFactory::release(preBasis);
      Traits::InterpolationFactory::release(interpol);

#if GLFEM_BASIS_PRINT
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,1> > BasisMI;
        const MBasisMI &_mBasisMI = Dune::MonomialBasisProvider<dimension,MIField>::template create<Topology>(order);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order << ".basis";
        std::ofstream out(name.str().c_str());
        basisPrint<0>(out,basisMI);
      }
#endif
      return basis;
    }
  };
}

#endif // #ifndef DUNE_LAGRANGEBASIS_HH
