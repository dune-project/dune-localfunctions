// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_L2INTERPOLATION_HH
#define DUNE_L2INTERPOLATION_HH

#include <dune/common/concept.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/utility/lfematrix.hh>

namespace Dune
{
  /**
   * @brief A local L2 interpolation taking a test basis and a quadrature
   *        rule.
   *
   * This class computes a local interpolation where the coefficients
   * are of the form:
   *     c = M^{-1}b
   * - M is the mass matrix with respect to the given basis and
   * - b = int f phi (where phi are the basis functions).
   * Thus the resulting local function u=c.varphi is defined through
   * the l2 interpolation int u phi = in f phi for all phi in the
   * base function set.
   * The third template argument can be used to specify that the
   * mass matrix is the unit matrix (onb=true).
   **/
  template< class B, class Q, bool onb >
  struct LocalL2Interpolation;

  template< class B, class Q >
  class LocalL2InterpolationBase
  {
    typedef LocalL2InterpolationBase< B, Q > This;

  public:
    typedef B Basis;
    typedef Q Quadrature;

    static const unsigned int dimension = Basis::dimension;

    /** \brief Interpolate a function that implements void evaluate(Domain, Range&) */
    template< class Function, class DofField, std::enable_if_t<models<Impl::FunctionWithEvaluate<typename Function::DomainType, typename Function::RangeType>, Function>(), int> = 0 >
    void interpolate ( const Function &function, std::vector< DofField > &coefficients ) const
    {
      typedef typename Quadrature::iterator Iterator;
      typedef FieldVector< DofField, Basis::dimRange > RangeVector;

      const unsigned int size = basis().size();
      static std::vector< RangeVector > basisValues( size );

      coefficients.resize( size );
      basisValues.resize( size );
      for( unsigned int i = 0; i < size; ++i )
        coefficients[ i ] = Zero< DofField >();

      const Iterator end = quadrature().end();
      for( Iterator it = quadrature().begin(); it != end; ++it )
      {
        basis().evaluate( it->position(), basisValues );
        typename Function::RangeType val;
        function.evaluate( field_cast<typename Function::DomainType::field_type>(it->position()), val );
        RangeVector factor = field_cast< DofField >( val );
        factor *= field_cast< DofField >( it->weight() );
        for( unsigned int i = 0; i < size; ++i )
          coefficients[ i ] += factor * basisValues[ i ];
      }
    }

    /** \brief Interpolate a function that implements Range operator()(Domain) */
    template< class Function, class DofField, std::enable_if_t<models<Impl::FunctionWithCallOperator<typename Quadrature::value_type::Vector>, Function>(), int> = 0 >
    void interpolate ( const Function &function, std::vector< DofField > &coefficients ) const
    {
      typedef FieldVector< DofField, Basis::dimRange > RangeVector;

      const unsigned int size = basis().size();
      static std::vector< RangeVector > basisValues( size );

      coefficients.resize( size );
      basisValues.resize( size );
      for( unsigned int i = 0; i < size; ++i )
        coefficients[ i ] = Zero< DofField >();

      for (auto&& qp : quadrature())
      {
        basis().evaluate( qp.position(), basisValues );
        auto val = function( qp.position() );
        RangeVector factor = field_cast< DofField >( val );
        factor *= field_cast< DofField >( qp.weight() );
        for( unsigned int i = 0; i < size; ++i )
          coefficients[ i ] += factor * basisValues[ i ];
      }
    }

    const Basis &basis () const
    {
      return basis_;
    }

    const Quadrature &quadrature () const
    {
      return quadrature_;
    }

  protected:
    LocalL2InterpolationBase ( const Basis &basis, const Quadrature &quadrature )
      : basis_( basis ),
        quadrature_( quadrature )
    {}

    const Basis &basis_;
    const Quadrature &quadrature_;
  };

  template< class B, class Q >
  struct LocalL2Interpolation<B,Q,true>
    : public LocalL2InterpolationBase<B,Q>
  {
    typedef LocalL2InterpolationBase<B,Q> Base;
    template< class BasisFactory, bool onb >
    friend class LocalL2InterpolationFactory;
    using typename Base::Basis;
    using typename Base::Quadrature;
  private:
    LocalL2Interpolation ( const typename Base::Basis &basis, const typename Base::Quadrature &quadrature )
      : Base(basis,quadrature)
    {}
  };
  template< class B, class Q >
  struct LocalL2Interpolation<B,Q,false>
    : public LocalL2InterpolationBase<B,Q>
  {
    typedef LocalL2InterpolationBase<B,Q> Base;
    template< class BasisFactory, bool onb >
    friend class LocalL2InterpolationFactory;
    using typename Base::Basis;
    using typename Base::Quadrature;
    template< class Function, class DofField >
    void interpolate ( const Function &function, std::vector< DofField > &coefficients ) const
    {
      const unsigned size = Base::basis().size();
      Base::interpolate(function,val_);
      coefficients.resize( size );
      for (unsigned int i=0; i<size; ++i)
      {
        coefficients[i] = 0;
        for (unsigned int j=0; j<size; ++j)
        {
          coefficients[i] += field_cast<DofField>(massMatrix_(i,j)*val_[j]);
        }
      }
    }
  private:
    LocalL2Interpolation ( const typename Base::Basis &basis, const typename Base::Quadrature &quadrature )
      : Base(basis,quadrature),
        val_(basis.size()),
        massMatrix_()
    {
      typedef FieldVector< Field, Base::Basis::dimRange > RangeVector;
      typedef typename Base::Quadrature::iterator Iterator;
      const unsigned size = basis.size();
      std::vector< RangeVector > basisValues( size );

      massMatrix_.resize( size,size );
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int j=0; j<size; ++j)
          massMatrix_(i,j) = 0;
      const Iterator end = Base::quadrature().end();
      for( Iterator it = Base::quadrature().begin(); it != end; ++it )
      {
        Base::basis().evaluate( it->position(), basisValues );
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size; ++j)
            massMatrix_(i,j) += (basisValues[i]*basisValues[j])*it->weight();
      }
      if ( !massMatrix_.invert() )
      {
        DUNE_THROW(MathError, "Mass matrix singular in LocalL2Interpolation");
      }

    }
    typedef typename Base::Basis::StorageField Field;
    typedef FieldVector< Field, Base::Basis::dimRange > RangeVector;
    typedef LFEMatrix<Field> MassMatrix;
    mutable std::vector<Field> val_;
    MassMatrix massMatrix_;
  };

  /**
   * @brief A factory class for the local l2 interpolations
   *        taking a basis factory.
   **/
  template< class BasisFactory, bool onb >
  struct LocalL2InterpolationFactory
  {
    static const unsigned int dimension = BasisFactory::dimension;
    typedef typename BasisFactory::Key Key;
    typedef typename BasisFactory::Object Basis;
    typedef double Field;
    typedef QuadratureRule<Field,dimension> Quadrature;
    typedef QuadratureRules<Field,dimension> QuadratureProvider;
    typedef LocalL2Interpolation< Basis, Quadrature, onb > LocalInterpolation;
    typedef const LocalInterpolation Object;

    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      constexpr Dune::GeometryType geometry = geometryId;
      const Basis *basis = BasisFactory::template create< geometry >( key );
      const Quadrature & quadrature = QuadratureProvider::rule(geometry, 2*basis->order()+1);
      return new Object( *basis, quadrature );
    }
    static void release ( Object *object )
    {
      const Basis &basis = object->basis();
      BasisFactory::release( &basis );
      delete object;
    }
  };

}

#endif // #ifndef DUNE_L2INTERPOLATION_HH
