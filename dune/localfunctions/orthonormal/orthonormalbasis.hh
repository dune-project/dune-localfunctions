// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH

#include <sstream>

#include <dune/localfunctions/utility/polynomialbasis.hh>
#include <dune/localfunctions/orthonormal/orthonormalcompute.hh>

namespace Dune
{

  // OrthonormalBasisFactory
  // -----------------------
  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisFactory
  {
    static const unsigned int dimension = dim;
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
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;

    typedef unsigned int Key;
    typedef const Basis Object;

    static constexpr GeometryType SimplexGeometry = GeometryTypes::simplex(dim);

    template< GeometryType::Id geometryId >
    static Object *create ( const unsigned int order )
    {
      const MonomialBasisType &monomialBasis = *MonomialBasisProviderType::template create< SimplexGeometry >( order );

      static CoefficientMatrix _coeffs;
      if( _coeffs.size() <= monomialBasis.size() )
      {
        ONBCompute::ONBMatrix< geometryId, ComputeField > matrix( order );
        _coeffs.fill( matrix );
      }

      return new Basis( monomialBasis, _coeffs, monomialBasis.size() );
    }
    static void release( Object *object ) { delete object; }
  };

}

#endif // #ifndef DUNE_ORTHONORMALBASIS_HH
