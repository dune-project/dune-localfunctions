// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEFAULTBASISFACTORY_HH
#define DUNE_DEFAULTBASISFACTORY_HH

#include <fstream>
#include <dune/common/exceptions.hh>
#include <dune/geometry/topologyfactory.hh>

#include <dune/localfunctions/utility/basismatrix.hh>

namespace Dune
{
  struct Identity
  {
    template <class T>
    static T apply( const T &t )
    {
      return t;
    }
  };
  /************************************************
  * Class for providing a factory for basis
  * functions over the set of reference elements.
  * Is based on the TopologyFactory but additionally
  * provides rebindes of the field type.
  * The user provides factories for the pre basis and the
  * interpolations. The default construction process of
  * the basis is performed in this class.
  ************************************************/
  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF,
      class PreBasisKeyExtractor = Identity >
  struct DefaultBasisFactory;

  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF,
      class PreBasisKeyExtractor >
  struct DefaultBasisFactoryTraits
  {
    static const unsigned int dimension = dim;
    static const unsigned int dimRange  = dimR;

    typedef PreBFactory PreBasisFactory;
    typedef typename PreBasisFactory::Object PreBasis;
    typedef InterpolFactory InterpolationFactory;
    typedef typename InterpolationFactory::Object Interpolation;

    typedef typename PreBasisFactory::template EvaluationBasisFactory<dim,SF>::Type MonomialBasisFactory;
    typedef typename MonomialBasisFactory::Object MonomialBasis;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< SF, dimRange > > Basis;

    typedef const Basis Object;
    typedef typename InterpolationFactory::Key Key;
    typedef DefaultBasisFactory<PreBFactory,InterpolFactory,dim,dimR,SF,CF,PreBasisKeyExtractor> Factory;
  };

  template< class PreBFactory,
      class InterpolFactory,
      unsigned int dim, unsigned int dimR,
      class SF, class CF,
      class PreBasisKeyExtractor >
  struct DefaultBasisFactory
    : public TopologyFactory<
          DefaultBasisFactoryTraits< PreBFactory,InterpolFactory,dim,dimR,SF,CF,PreBasisKeyExtractor >
          >
  {
    typedef DefaultBasisFactoryTraits< PreBFactory,InterpolFactory,dim,dimR,SF,CF,PreBasisKeyExtractor > Traits;
    static const unsigned int dimension = Traits::dimension;
    static const unsigned int dimRange  = Traits::dimRange;
    typedef SF StorageField;
    typedef CF ComputeField;
    typedef typename Traits::Basis Basis;
    typedef typename Traits::PreBasisFactory PreBasisFactory;

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
      const typename PreBasisFactory::Key preBasisKey = PreBasisKeyExtractor::apply(key);
      const typename Traits::PreBasis *preBasis = Traits::PreBasisFactory::template create<Topology>( preBasisKey );
      const typename Traits::Interpolation *interpol = Traits::InterpolationFactory::template create<Topology>( key );
      BasisMatrix< typename Traits::PreBasis,
          typename Traits::Interpolation,
          ComputeField > matrix( *preBasis, *interpol );

      const typename Traits::MonomialBasis *monomialBasis = Traits::MonomialBasisFactory::template create< Topology >( preBasis->order() );

      Basis *basis = new Basis( *monomialBasis );

      basis->fill( matrix );

      Traits::InterpolationFactory::release(interpol);
      Traits::PreBasisFactory::release(preBasis);

      return basis;
    }
    //! release the object returned by the create methods
    static void release( Object *object)
    {
      const typename Traits::MonomialBasis *monomialBasis = &(object->basis());
      delete object;
      Traits::MonomialBasisFactory::release( monomialBasis );
    }
  };
}

#endif // #ifndef DUNE_DEFAULTBASISFACTORY_HH
