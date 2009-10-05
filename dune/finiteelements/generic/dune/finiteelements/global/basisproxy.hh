// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FINITEELEMENTS_BASISPROXY_HH
#define DUNE_FINITEELEMENTS_BASISPROXY_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  template< class Basis, class Geometry >
  class BasisProxy
  {
    typedef BasisProxy< Basis, Geometry > This;

  public:
    typedef typename Basis::Field Field;
    static const int dimRange = Basis::dimRange;
    static const int dimDomain = Basis::dimension;
    static const int dimWorld = Geometry::coorddimension;
    typedef typename Basis::DomainVector DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;
    typedef FieldMatrix< Field, dimRange, dimWorld > Jacobian;

    BasisProxy ( const Basis &basis, const Geometry &geometry )
      : basis_( &basis ),
        geometry_( &geometry )
    {}

    const Basis &basis () const
    {
      return *basis_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    const unsigned int order () const
    {
      return basis().order();
    }

    const unsigned int size () const
    {
      return basis().size();
    }

    template <class Fx>
    void evaluate ( const FieldVector<Fx,dimDomain> &x, std::vector< RangeVector > &values ) const
    {
      return evaluate( field_cast<Field>(x),values );
    }
    template< class Fx, class LocalDofVector >
    void evaluate ( const FieldVector<Fx,dimDomain>&x, const LocalDofVector &localDofs, RangeVector &value ) const
    {
      return evaluate( field_cast<Field>(x),localDofs,value );
    }

    void evaluate ( const DomainVector &x, std::vector< RangeVector > &values ) const
    {
      values.resize( size() );
      return basis().evaluate( x, values );
    }

    template< class LocalDofVector >
    void evaluate ( const DomainVector &x, const LocalDofVector &localDofs, RangeVector &value ) const
    {
      static std::vector< RangeVector > values;
      evaluate( x, values );

      value = Zero< Field >();
      for( unsigned int i = 0; i < size(); ++i )
        value.axpy( localDofs[ i ], values[ i ] );
    }

    void jacobian ( const DomainVector &x, std::vector< Jacobian > &jacobians ) const
    {
      typedef Dune::Derivatives< Field, dimDomain, dimRange, 1, value > Derivatives;
      static std::vector< Derivatives > derivatives;
      derivatives.resize( basis.size() );
      basis().template evaluate< 1 >( x, derivatives );
      const Dune::FieldMatrix< Field, dimWorld, dimDomain > gjit = geometry().jacobianInterseTransposed( x );
      for( unsigned int i = 0; i < basis.size(); ++i )
      {
        for( unsigned int r = 0; r < dimRange; ++r )
          gjit.mv( derivatives[ i ][ r ], jacobians[ i ][ r ] );
      }
    }

    template< class LocalDofVector >
    void jacobian ( const DomainVector &x, const LocalDofVector &localDofs, Jacobian &jacobian ) const
    {
      typedef Dune::Derivatives< Field, dimDomain, dimRange, 1, value > Derivatives;
      static std::vector< Derivatives > derivatives;
      derivatives.resize( basis.size() );
      basis().template evaluate< 1 >( x, derivatives );
      Dune::FieldMatrix< Field, dimRange, dimDomain > j( Zero< Field >() );
      for( unsigned int i = 0; i < basis.size(); ++i )
      {
        for( unsigned int r = 0; r < dimRange; ++r )
          j[ r ].axpy( localDofs[ i ], derivatives[ i ][ r ] );
      }
      const Dune::FieldMatrix< Field, dimWorld, dimDomain > gjit = geometry().jacobianInterseTransposed( x );
      for( unsigned int r = 0; r < dimRange; ++r )
        gjit.mv( j[ r ], jacobian[ r ] );
    }

  private:
    const Basis *basis_;
    const Geometry *geometry_;
  };

}

#endif // #ifndef DUNE_FINITEELEMENTS_BASISPROXY_HH
