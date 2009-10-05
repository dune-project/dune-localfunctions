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
    typedef typename Basis::DomainVector DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;
    typedef FieldMatrix< Field, dimRange, dimDomain > Jacobian;

    BasisProxy ( const Basis &basis, const Geometry &geoetry )
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
      // ...
    }

    template< class LocalDofVector >
    void jacobian ( const DomainVector &x, const LocalDofVector &localDofs, Jacobian &jacobian ) const
    {
      // ...
    }

  private:
    const Basis *basis_;
    const Geometry *geometry_;
  };

}

#endif // #ifndef DUNE_FINITEELEMENTS_BASISPROXY_HH
