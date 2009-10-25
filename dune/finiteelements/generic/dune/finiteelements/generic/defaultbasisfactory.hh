// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEFAULTBASISFACTORY_HH
#define DUNE_DEFAULTBASISFACTORY_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/finiteelements/generic/basismatrix.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>

namespace Dune
{
  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF >
  struct DefaultBasisFactory;

  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF >
  struct DefaultBasisFactoryTraits
  {
    static const unsigned int dimension = dim;
    static const unsigned int dimRange  = dimR;

    typedef PreBFactory PreBasisFactory;
    typedef typename PreBasisFactory::Object PreBasis;
    typedef InterpolFactory InterpolationFactory;
    typedef typename InterpolationFactory::Object Interpolation;

    typedef typename PreBasisFactory::template EvaluationBasisFactory<dim,SF>::Type MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Object MonomialBasis;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< SF, dimRange > > Basis;

    typedef const Basis Object;
    typedef typename PreBasisFactory::Key Key;  // should be more flexible
    typedef DefaultBasisFactory<PreBFactory,InterpolFactory,dim,dimR,SF,CF> Factory;
  };

  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF >
  struct DefaultBasisFactory
    : public TopologyFactory<
          DefaultBasisFactoryTraits< PreBFactory,InterpolFactory,dim,dimR,SF,CF >
          >
  {
    typedef DefaultBasisFactoryTraits< PreBFactory,InterpolFactory,dim,dimR,SF,CF > Traits;
    static const unsigned int dimension = dim;
    typedef SF StorageField;
    typedef CF ComputeField;
    typedef typename Traits::Basis Basis;

    typedef typename Traits::Object Object;
    typedef typename Traits::Key Key;
    template <unsigned int dd, class FF>
    struct EvaluationBasisFactory
    {
      typedef typename Traits::PreBasisFactory::template EvaluationBasisFactory<dd,FF>::Type
      Type;
    };

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      const typename Traits::PreBasis *preBasis = Traits::PreBasisFactory::template create<Topology>( key );
      const typename Traits::Interpolation *interpol = Traits::InterpolationFactory::template create<Topology>( key );
      BasisMatrix< typename Traits::PreBasis,
          typename Traits::Interpolation,
          ComputeField > matrix( *preBasis, *interpol );

      const typename Traits::MonomialBasis *monomialBasis = Traits::MonomialBasisProvider::template create< Topology >( preBasis->order() );

      Basis *basis = new Basis( *monomialBasis );

      basis->fill( matrix );

      Traits::PreBasisFactory::release(preBasis);
      Traits::InterpolationFactory::release(interpol);

      return basis;
    }
  };
}

#endif // #ifndef DUNE_DEFAULTBASISFACTORY_HH
