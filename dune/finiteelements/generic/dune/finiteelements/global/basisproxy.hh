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
    typedef typename Basis::DomainVector DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

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
      return basis().evaluate( x, values );
    }

  private:
    const Basis *basis_;
    const Geometry *geometry_;
  };

}

#endif // #ifndef DUNE_FINITEELEMENTS_BASISPROXY_HH
