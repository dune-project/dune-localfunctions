// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/misc.hh>

#include <dune/common/field.hh>

#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/multiindex.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology >
  class MonomialBasisSize;

  template< class Topology, class F >
  class MonomialBasis;



  // MonomialBasisSize
  // -----------------

  template<>
  class MonomialBasisSize< GenericGeometry::Point >
  {
  public:
    typedef MonomialBasisSize< GenericGeometry::Point > This;
    static This &instance()
    {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Point Topology;
    static const unsigned int dimDomain = Topology::dimension;

    friend class MonomialBasisSize< GenericGeometry::Prism< Topology > >;
    friend class MonomialBasisSize< GenericGeometry::Pyramid< Topology > >;

    mutable unsigned int maxOrder_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisSize ()
      : maxOrder_( 0 ),
        sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisSize ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    unsigned int maxOrder () const
    {
      return maxOrder_;
    }

    void computeSizes ( unsigned int order ) const
    {
      if (order <= maxOrder_)
        return;

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

  template< class BaseTopology >
  class MonomialBasisSize< GenericGeometry::Prism< BaseTopology > >
  {

  public:
    typedef MonomialBasisSize< GenericGeometry::Prism< BaseTopology > > This;
    static This &instance() {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Prism< BaseTopology > Topology;
    static const unsigned int dimDomain = Topology::dimension;

    friend class MonomialBasisSize< GenericGeometry::Prism< Topology > >;
    friend class MonomialBasisSize< GenericGeometry::Pyramid< Topology > >;

    mutable unsigned int maxOrder_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisSize ()
      : maxOrder_( 0 ),
        sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisSize ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    unsigned int maxOrder() const
    {
      return maxOrder_;
    }

    void computeSizes ( unsigned int order ) const
    {
      if (order <= maxOrder_)
        return;

      maxOrder_ = order;

      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      MonomialBasisSize<BaseTopology> &baseBasis =
        MonomialBasisSize<BaseTopology>::instance();
      baseBasis.computeSizes( order );
      const unsigned int *const baseSizes = baseBasis.sizes_;
      const unsigned int *const baseNBF   = baseBasis.numBaseFunctions_;

      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ] + k*baseSizes[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
      }
    }
  };

  template< class BaseTopology >
  class MonomialBasisSize< GenericGeometry::Pyramid< BaseTopology > >
  {
  public:
    typedef MonomialBasisSize< GenericGeometry::Pyramid< BaseTopology > > This;
    static This &instance() {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Pyramid< BaseTopology > Topology;
    static const unsigned int dimDomain = Topology::dimension;

    friend class MonomialBasisSize< GenericGeometry::Prism< Topology > >;
    friend class MonomialBasisSize< GenericGeometry::Pyramid< Topology > >;

    mutable unsigned int maxOrder_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisSize ()
      : maxOrder_( 0 ),
        sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisSize ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    unsigned int maxOrder() const
    {
      return maxOrder_;
    }

    void computeSizes ( unsigned int order ) const
    {
      if (order <= maxOrder_)
        return;

      maxOrder_ = order;

      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      MonomialBasisSize<BaseTopology> &baseBasis =
        MonomialBasisSize<BaseTopology>::instance();

      baseBasis.computeSizes( order );

      const unsigned int *const baseNBF = baseBasis.numBaseFunctions_;
      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
      }
    }
  };

  // Internal Forward Declarations
  // -----------------------------


  template< int dim, class F >
  struct MonomialBasisHelper
  {
    typedef MonomialBasisSize< typename GenericGeometry::SimplexTopology< dim >::type > Size;

    static void
    copy ( const unsigned int deriv, F *&wit, F *&rit, const unsigned int numBaseFunc, const F &z )
    {
      Size &size = Size::instance();
      const F *const rend = rit + numBaseFunc;
      for( ; rit != rend; )
      {
        for ( unsigned d=0; d<=deriv; ++d )
        {
          const unsigned int endS = size.sizes_[d];
          for ( unsigned s=0; s<endS; ++s,++rit, ++wit )
            *wit = z * *rit;
        }
      }
    }
    static void
    setPoint( const unsigned int deriv, F *&wit)
    {
      // Size &size = Size::instance();
      *wit = Unity<F>();
    }
  };



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

    typedef MonomialBasisHelper< dimDomain, Field > Helper;

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      Field *row0 = values;
      Helper::setPoint(deriv,row0);
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      values[ 0 ] = Unity< Field >();
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

    typedef MonomialBasisHelper< dimDomain, Field > Helper;
    typedef MonomialBasisSize< BaseTopology > BaseSize;
    typedef MonomialBasisSize< Topology > Size;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;

    MonomialBasisImpl ()
    {}

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();

      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( deriv , order, x, offsets, values );

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *row1 = values + offsets[ k-1 ];
        Field *wit = row1 + size.sizes_[ k ];
        Helper::copy( deriv, wit, row1, k*size.sizes_[ k ], z );
        Helper::copy( deriv, wit, row0, size.numBaseFunctions_[ k-1 ], z );
        row0 = row1;
      }
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();
      const Size &mySize = Size::instance();
      // fill first column
      baseBasis_.integral( order, offsets, values );
      const unsigned int *const baseSizes = size.sizes_;

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1begin = values + offsets[ k-1 ];
        Field *const row1End = row1begin + mySize.sizes_[ k ];
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

    typedef MonomialBasisHelper< dimDomain, Field > Helper;
    typedef MonomialBasisSize< BaseTopology > BaseSize;
    typedef MonomialBasisSize< Topology > Size;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;

    MonomialBasisImpl ()
    {}

    template< int dimD >
    void evaluateSimplexBase ( const unsigned int deriv, const unsigned int order,
                               const FieldVector< Field, dimD > &x,
                               const unsigned int *const offsets,
                               Field *const values,
                               const BaseSize &size ) const
    {
      baseBasis_.evaluate( deriv, order, x, offsets, values );
    }

    template< int dimD >
    void evaluatePyramidBase ( const unsigned int deriv, const unsigned int order,
                               const FieldVector< Field, dimD > &x,
                               const unsigned int *const offsets,
                               Field *const values,
                               const BaseSize &size ) const
    {

      Field omz = Unity< Field >() - x[ dimDomain-1 ];

      if( Zero< Field >() < omz )
      {
        const Field invomz = Unity< Field >() / omz;
        FieldVector< Field, dimD > y;
        for( unsigned int i = 0; i < dimDomain-1; ++i )
          y[ i ] = x[ i ] * invomz;

        // fill first column
        baseBasis_.evaluate( deriv, order, y, offsets, values );

        Field omzk = omz;
        for( unsigned int k = 1; k <= order; ++k )
        {
          Field *it = values + offsets[ k-1 ];
          Field *const end = it + size.sizes_[ k ];
          for( ; it != end; ++it )
            *it *= omzk;
          omzk *= omz;
        }
      }
      else
      {
        assert( deriv==0 );
        *values = Unity< Field >();
        for( unsigned int k = 1; k <= order; ++k )
        {
          Field *it = values + offsets[ k-1 ];
          Field *const end = it + size.sizes_[ k ];
          for( ; it != end; ++it )
            *it = Zero< Field >();
        }
      }
    }

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();

      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplexBase( deriv, order, x, offsets, values, size );
      else
        evaluatePyramidBase( deriv, order, x, offsets, values, size );

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *row1 = values + offsets[ k-1 ];
        Field *wit = row1 + size.sizes_[ k ];
        Helper::copy( deriv, wit, row0, size.numBaseFunctions_[ k-1 ], x[ dimDomain-1 ] );
        row0 = row1;
      }
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();

      // fill first column
      baseBasis_.integral( order, offsets, values );

      const unsigned int *const baseSizes = size.sizes_;

      Field *const col0End = values + baseSizes[ 0 ];
      for( Field *it = values; it != col0End; ++it )
        *it *= Field( 1 ) /  Field( double(dimDomain) );  // ??? double cast due to error in Linker
      Field *row0 = values;

      for( unsigned int k = 1; k <= order; ++k )
      {
        const Field factor = (Field( 1 ) / Field( k + dimDomain ));

        Field *const row1 = values+offsets[ k-1 ];
        Field *const col0End = row1 + baseSizes[ k ];
        Field *it = row1;
        for( ; it != col0End; ++it )
          *it *= factor;
        for( unsigned int i = 1; i <= k; ++i )
        {
          Field *const end = it + baseSizes[ k-i ];
          assert( (unsigned int)(end - values) <= offsets[ k ] );
          for( ; it != end; ++row0, ++it )
            *it = (*row0) * (Field( i ) * factor);
        }
        row0 = row1;
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
    static const unsigned int dimension = Base::dimDomain;

    typedef typename Base::Field Field;

    typedef typename Base::DomainVector DomainVector;

    typedef MonomialBasisSize<Topology> Size;

    MonomialBasis (unsigned int order)
      : Base(),
        order_(order),
        size_(Size::instance())
    {}

    const unsigned int *sizes ( unsigned int order ) const
    {
      size_.computeSizes( order );
      return size_.numBaseFunctions_;
    }
    const unsigned int size ( ) const
    {
      return sizes( order_ )[ order_ ];
    }
    const unsigned int order ( ) const
    {
      return order_ ;
    }

    void evaluate ( unsigned int deriv, const DomainVector &x,
                    Field *const values ) const
    {
      typedef typename GenericGeometry::SimplexTopology< dimension >::type SimplexTopology;
      MonomialBasisSize< SimplexTopology >::instance().computeSizes( deriv );
      Base::evaluate( deriv, order_, x, sizes( order_ ), values );
    }

    template <unsigned int deriv>
    void evaluate ( const DomainVector &x,
                    Field *const values ) const
    {
      evaluate( deriv, x, values );
    }

    template <unsigned int deriv, class F1 >
    void evaluate ( const DomainVector &x,
                    F1 *const values ) const
    {
      evaluate<deriv>( x, reinterpret_cast< Field * >( values ) );
    }
    template<unsigned int deriv, class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<deriv>(x,&(values[0]));
    }
    template <class F1>
    void evaluate ( const DomainVector &x,
                    F1 *const values ) const
    {
      evaluate<0>( x, reinterpret_cast< Field * >( values ) );
    }
    template<class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<0>(x,&(values[0]));
    }

    void integral ( Field *const values ) const
    {
      Base::integral( order_, sizes( order_ ), values );
    }
    template <class RangeVector>
    void integral ( std::vector< RangeVector > &values ) const
    {
      integral( &(values[ 0 ]) );
    }
  private:
    MonomialBasis(const This&);
    This& operator=(const This&);
    unsigned int order_;
    Size &size_;
  };



  // StdMonomialBasis
  // ----------------

  template< int dim,class F >
  class StandardMonomialBasis
    : public MonomialBasis< typename GenericGeometry::SimplexTopology< dim >::type, F >
  {
    typedef StandardMonomialBasis< dim, F > This;
    typedef MonomialBasis< typename GenericGeometry::SimplexTopology< dim >::type, F > Base;

  public:
    typedef typename GenericGeometry::SimplexTopology< dim >::type Topology;
    static const int dimension = dim;

    StandardMonomialBasis ( unsigned int order )
      : Base( order )
    {}
  };



  // StandardBiMonomialBasis
  // -----------------------

  template< int dim, class F >
  class StandardBiMonomialBasis
    : public MonomialBasis< typename GenericGeometry::CubeTopology< dim >::type, F >
  {
    typedef StandardBiMonomialBasis< dim, F > This;
    typedef MonomialBasis< typename GenericGeometry::CubeTopology< dim >::type, F > Base;

  public:
    typedef typename GenericGeometry::CubeTopology< dim >::type Topology;
    static const int dimension = dim;

    StandardBiMonomialBasis ( unsigned int order )
      : Base( order )
    {}
  };

  // -----------------------------------------------------------
  // -----------------------------------------------------------
  // VirtualMonomialBasis
  // -------------------

  template< int dim, class F >
  class VirtualMonomialBasis // : public BasisInterface<dim,F>
  {
    typedef VirtualMonomialBasis< dim, F > This;

  public:
    typedef F Field;
    static const int dimension = dim;

    typedef FieldVector<Field,dimension> DomainVector;

    explicit VirtualMonomialBasis(unsigned int order)
      : order_(order) {}

    virtual ~VirtualMonomialBasis() {}

    virtual const unsigned int *sizes ( ) const = 0;

    const unsigned int size ( ) const
    {
      return sizes( )[ order_ ];
    }

    const unsigned int order () const
    {
      return order_;
    }

    virtual void evaluate ( const unsigned int deriv, const DomainVector &x,
                            Field *const values ) const = 0;
    template < unsigned int deriv, class Field,int size >
    void evaluate ( const DomainVector &x,
                    Dune::FieldVector<Field,size> *const values ) const
    {
      evaluate( deriv, x, reinterpret_cast< Field * >( values ) );
    }
    template <unsigned int deriv, class RangeVector>
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( deriv, x, &(values[ 0 ]) );
    }
    template< class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<0>(x,values);
    }

    virtual void integral ( Field *const values ) const = 0;
    template <class RangeVector>
    void integral ( std::vector< RangeVector > &values ) const
    {
      integral( &(values[ 0 ]) );
    }
  protected:
    unsigned int order_;
  };

  template< class Topology, class F >
  class VirtualMonomialBasisImpl
    : public VirtualMonomialBasis< Topology::dimension, F >
  {
    typedef VirtualMonomialBasis< Topology::dimension, F > Base;
    typedef VirtualMonomialBasisImpl< Topology, F > This;

  public:
    typedef typename Base::Field Field;
    typedef typename Base::DomainVector DomainVector;

    VirtualMonomialBasisImpl(unsigned int order)
      : Base(order), basis_(order)
    {}

    const unsigned int *sizes ( ) const
    {
      return basis_.sizes(order_);
    }

    void evaluate ( const unsigned int deriv, const DomainVector &x,
                    Field *const values ) const
    {
      basis_.evaluate(deriv,x,values);
    }

    void integral ( Field *const values ) const
    {
      basis_.integral(values);
    }

  private:
    MonomialBasis<Topology,Field> basis_;
    using Base::order_;
  };

  template< int dim, class F >
  struct MonomialBasisCreator
  {
    typedef F StorageField;
    typedef VirtualMonomialBasis<dim,StorageField> Basis;
    static const int dimension = dim;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        basis = new VirtualMonomialBasisImpl<Topology,StorageField>(order);
      }
    };
  };

  template< int dim, class SF >
  struct MonomialBasisProvider
    : public BasisProvider<MonomialBasisCreator<dim,SF> >
  {};

}

#endif
