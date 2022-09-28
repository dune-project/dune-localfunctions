// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DEFAULTBASISFACTORY_HH
#define DUNE_DEFAULTBASISFACTORY_HH

#include <fstream>
#include <dune/common/exceptions.hh>

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
  struct DefaultBasisFactory
  {
    static const unsigned int dimension = dim;
    static const unsigned int dimRange  = dimR;
    typedef SF StorageField;
    typedef CF ComputeField;
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
    template <unsigned int dd, class FF>
    struct EvaluationBasisFactory
    {
      typedef typename PreBasisFactory::template EvaluationBasisFactory<dd,FF>::Type
      Type;
    };

    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      const typename PreBasisFactory::Key preBasisKey = PreBasisKeyExtractor::apply(key);
      const PreBasis *preBasis = PreBasisFactory::template create<geometryId>( preBasisKey );
      const Interpolation *interpol = InterpolationFactory::template create<geometryId>( key );
      BasisMatrix< PreBasis, Interpolation, ComputeField > matrix( *preBasis, *interpol );

      const MonomialBasis *monomialBasis = MonomialBasisFactory::template create< geometryId >( preBasis->order() );

      Basis *basis = new Basis( *monomialBasis );

      basis->fill( matrix );

      InterpolationFactory::release(interpol);
      PreBasisFactory::release(preBasis);

      return basis;
    }
    //! release the object returned by the create methods
    static void release( Object *object)
    {
      const MonomialBasis *monomialBasis = &(object->basis());
      delete object;
      MonomialBasisFactory::release( monomialBasis );
    }
  };
}

#endif // #ifndef DUNE_DEFAULTBASISFACTORY_HH
