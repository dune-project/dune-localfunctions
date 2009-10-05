// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology, class F >
  class MonomialBasis;



  // MonomialBasisImpl
  // -----------------

  template< class Topology, class F >
  class MonomialBasisImpl;

  template< class F >
  class MonomialBasisImpl< GenericGeometry::Point, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Point, F > This;

  public:
    typedef GenericGeometry::Point Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    mutable unsigned int maxOrder_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      values[ 0 ] = Field( 1 );
    }

    unsigned int maxOrder () const
    {
      return maxOrder_;
    }

    void computeSizes ( unsigned int order ) const
    {
      maxOrder_ = order;

      delete [] sizes_;
      delete [] numBaseFunctions_;
      sizes_            = new unsigned int [ order+1 ];
      numBaseFunctions_ = new unsigned int [ order+1 ];

      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = 0;
        numBaseFunctions_[ k ] = 1;
      }
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasisImpl< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Prism< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Prism< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );
      const unsigned int *const baseSizes = baseBasis_.sizes_;

      RangeVector *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        RangeVector *row1 = values + offsets[ k-1 ];
        RangeVector *it = row1 + baseSizes[ k ];
        RangeVector *const colkEnd = row1 + (k+1)*baseSizes[ k ];
        for( ; it != colkEnd; ++row1, ++it )
          *it = z * (*row1);
        RangeVector *const row1End = row1 + sizes_[ k ];
        for( ; it!=row1End; ++row0,++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    unsigned int maxOrder() const
    {
      return baseBasis_.maxOrder();
    }

    void computeSizes ( unsigned int order ) const
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      baseBasis_.computeSizes( order );
      const unsigned int *const baseSizes = baseBasis_.sizes_;
      const unsigned int *const baseNBF   = baseBasis_.numBaseFunctions_;

      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ] + k*baseSizes[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + baseNBF[ k ];
      }
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasisImpl< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Pyramid< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Pyramid< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluateSimplex ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           RangeVector *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );

      const unsigned int *const baseSizes = baseBasis_.sizes_;
      RangeVector *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        RangeVector *const row1 = values+offsets[ k-1 ];
        RangeVector *const row1End = row1+sizes_[ k ];
        for( RangeVector *it = row1 + baseSizes[ k ]; it!=row1End; ++row0,++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    template< int dimD >
    void evaluatePyramid ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           RangeVector *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];
      Field omz = Field( 1 ) - z;

      if( omz > 1e-12 ) // this number must depend on Field
      {
        const Field invomz = Field( 1 ) / omz;
        FieldVector< Field, dimDomain-1 > y;
        for( unsigned int i = 0; i < dimDomain-1; ++i )
          y[ i ] = x[ i ] * invomz;

        // fill first column
        baseBasis_.evaluate( order, y, offsets, values );
      }
      else
        omz = Field( 0 );

      const unsigned int *const baseSizes = baseBasis_.sizes_;
      RangeVector *row0 = values;
      Field omzk = omz;
      for( unsigned int k = 1; k <= order; ++k )
      {
        RangeVector *const row1 = values + offsets[ k-1 ];
        RangeVector *const row1End = row1 + sizes_[ k ];
        RangeVector *const col0End = row1 + baseSizes[ k ];
        RangeVector *it = row1;
        for( ; it!=col0End; ++it )
          *it = (*it) * omzk;
        for( ; it!=row1End; ++row0,++it )
          *it = z * (*row0);
        row0 = row1;
        omzk *= omz;
      }
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplex( order, x, offsets, values );
      else
        evaluatePyramid( order, x, offsets, values );
    }

    unsigned int maxOrder() const
    {
      return baseBasis_.maxOrder();
    }

    void computeSizes ( unsigned int order ) const
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      baseBasis_.computeSizes( order );

      const unsigned int *const baseNBF = baseBasis_.numBaseFunctions_;
      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + baseNBF[ k ];
      }
    }
  };



  // MonomialBasis
  // -------------

  template< class Topology, class F >
  class MonomialBasis
    : public MonomialBasisImpl< Topology, F >
  {
    typedef MonomialBasis< Topology, F > This;
    typedef MonomialBasisImpl< Topology, F > Base;

  public:
    typedef typename Base::Field Field;

    typedef typename Base::DomainVector DomainVector;
    typedef typename Base::RangeVector RangeVector;

    MonomialBasis ()
      : Base()
    {}

    const unsigned int *sizes ( unsigned int order ) const
    {
      if( order > Base::maxOrder() )
        Base::computeSizes( order );
      return Base::numBaseFunctions_;
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      Base::evaluate( order, x, sizes( order ), values );
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( order, x, &(values[ 0 ]) );
    }
  };

}

#endif
