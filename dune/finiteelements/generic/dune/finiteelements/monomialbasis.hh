// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/finiteelements/field.hh>

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

    typedef FieldVector< Field, dimDomain > DomainVector;

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
                    Field *const values ) const
    {
      values[ 0 ] = Unity<Field>();
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
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

    typedef FieldVector< Field, dimDomain > DomainVector;

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
                    Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );
      const unsigned int *const baseSizes = baseBasis_.sizes_;

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1begin = values + offsets[ k-1 ];
        Field *const colkEnd = row1begin + (k+1)*baseSizes[ k ];
        assert( (unsigned int)(colkEnd - values) <= offsets[ k ] );
        Field *const row1End = row1begin + sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );

        Field *row1 = row1begin;
        Field *it;
        for( it = row1begin + baseSizes[ k ]; it != colkEnd; ++row1, ++it )
          *it = z * (*row1);
        for( ; it != row1End; ++row0, ++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      // fill first column
      baseBasis_.integral( order, offsets, values );
      const unsigned int *const baseSizes = baseBasis_.sizes_;

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1begin = values + offsets[ k-1 ];
        Field *const row1End = row1begin + sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );

        Field *row1 = row1begin;
        Field *it = row1begin + baseSizes[ k ];
        for( unsigned int j = 1; j <= k; ++j )
        {
          Field *const end = it + baseSizes[ k ];
          assert( (unsigned int)(end - values) <= offsets[ k ] );
          for( ; it != end; ++row1, ++it )
            *it = (Field( j ) / Field( j+1 )) * (*row1);
        }
        for( ; it != row1End; ++row0, ++it )
          *it = (Field( k ) / Field( k+1 )) * (*row0);
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
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
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

    typedef FieldVector< Field, dimDomain > DomainVector;

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
                           Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );

      const unsigned int *const baseSizes = baseBasis_.sizes_;
      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1 = values+offsets[ k-1 ];
        Field *const row1End = row1+sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );
        for( Field *it = row1 + baseSizes[ k ]; it != row1End; ++row0, ++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    template< int dimD >
    void evaluatePyramid ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];
      Field omz = Unity< Field >() - z;

      if( Zero<Field>() < omz )
      {
        const Field invomz = Unity<Field>() / omz;
        FieldVector< Field, dimDomain-1 > y;
        for( unsigned int i = 0; i < dimDomain-1; ++i )
          y[ i ] = x[ i ] * invomz;

        // fill first column
        baseBasis_.evaluate( order, y, offsets, values );
      }
      else
        omz = Zero<Field>();

      const unsigned int *const baseSizes = baseBasis_.sizes_;
      Field *row0 = values;
      Field omzk = omz;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1 = values + offsets[ k-1 ];
        Field *const row1End = row1 + sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );
        Field *const col0End = row1 + baseSizes[ k ];
        Field *it = row1;
        for( ; it != col0End; ++it )
          *it = (*it) * omzk;
        for( ; it != row1End; ++row0, ++it )
          *it = z * (*row0);
        row0 = row1;
        omzk *= omz;
      }
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplex( order, x, offsets, values );
      else
        evaluatePyramid( order, x, offsets, values );
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      // fill first column
      baseBasis_.integral( order, offsets, values );

      const unsigned int *const baseSizes = baseBasis_.sizes_;

      Field *const col0End = values + baseSizes[ 0 ];
      for( Field *it = values; it != col0End; ++it )
        *it *= Field( 1 ) / Field( dimDomain );
      Field *row0 = values;

      for( unsigned int k = 1; k <= order; ++k )
      {
        const Field factor = (Field( 1 ) / Field( k + dimDomain ));

        Field *const row1 = values+offsets[ k-1 ];
        Field *const col0End = row1 + baseSizes[ k ];
        Field *it = row1;
        for( ; it != col0End; ++it )
          *it *= factor;
        for( unsigned int i = 0; i < k; ++i )
        {
          Field *const end = it + baseSizes[ i ];
          assert( (unsigned int)(end - values) <= offsets[ k ] );
          for( ; it != end; ++row0, ++it )
            *it = (factor * Field( i+1 )) * (*row0);
        }
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

      const unsigned int *const baseNBF = baseBasis_.numBaseFunctions_;
      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
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

    static const unsigned int dimRange = 1;

    typedef typename Base::DomainVector DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

    MonomialBasis ()
      : Base()
    {}

    const unsigned int *sizes ( unsigned int order ) const
    {
      if( order > Base::maxOrder() )
        Base::computeSizes( order );
      return Base::numBaseFunctions_;
    }

    const unsigned int size ( unsigned int order ) const
    {
      return sizes( order )[ order ];
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    Field *const values ) const
    {
      Base::evaluate( order, x, sizes( order ), values );
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      evaluate( order, x, reinterpret_cast< Field * >( values ) );
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( order, x, &(values[ 0 ]) );
    }

    void integral ( const unsigned int order,
                    Field *const values ) const
    {
      Base::integral( order, sizes( order ), values );
    }

    void integral ( const unsigned int order,
                    RangeVector *const values ) const
    {
      integral( order, reinterpret_cast< Field * >( values ) );
    }

    void integral ( const unsigned int order,
                    std::vector< RangeVector > &values ) const
    {
      integral( order, &(values[ 0 ]) );
    }

  };



  // StdMonomialTopology
  // -------------------

  template <int dim>
  struct StdMonomialTopology {
    typedef StdMonomialTopology<dim-1> BaseType;
    typedef GenericGeometry::Pyramid<typename BaseType::Type> Type;
  };
  template <>
  struct StdMonomialTopology<0> {
    typedef GenericGeometry::Point Type;
  };
  template <int dim>
  struct StdBiMonomialTopology {
    typedef GenericGeometry::Prism<typename StdBiMonomialTopology<dim-1>::Type> Type;
  };
  template <>
  struct StdBiMonomialTopology<0> {
    typedef GenericGeometry::Point Type;
  };
  template< int dim,class F >
  class StandardMonomialBasis
    : public MonomialBasis< typename StdMonomialTopology<dim>::Type, F >
  {
  public:
    typedef typename StdMonomialTopology<dim>::Type Topology;
  private:
    typedef StandardMonomialBasis< dim, F > This;
    typedef MonomialBasis< Topology, F > Base;
  public:
    StandardMonomialBasis ()
      : Base()
    {}
  };
  template< int dim, class F >
  class StandardBiMonomialBasis
    : public MonomialBasis< typename StdBiMonomialTopology<dim>::Type, F >
  {
  public:
    typedef typename StdBiMonomialTopology<dim>::Type Topology;
  private:
    typedef StandardBiMonomialBasis< dim, F > This;
    typedef MonomialBasis< Topology, F > Base;
  public:
    StandardBiMonomialBasis ()
      : Base()
    {}
  };

}

#endif
