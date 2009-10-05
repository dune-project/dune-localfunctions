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
#include <dune/finiteelements/tensor.hh>

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
    typedef MonomialBasisSize< GenericGeometry::Point > This;

  public:
    static This &instance ()
    {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Point Topology;

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

    unsigned int operator() ( const unsigned int order ) const
    {
      return numBaseFunctions_[ order ];
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
    typedef MonomialBasisSize< GenericGeometry::Prism< BaseTopology > > This;

  public:
    static This &instance ()
    {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Prism< BaseTopology > Topology;

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

    unsigned int operator() ( const unsigned int order ) const
    {
      return numBaseFunctions_[ order ];
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
    typedef MonomialBasisSize< GenericGeometry::Pyramid< BaseTopology > > This;

  public:
    static This &instance ()
    {
      static This _instance;
      return _instance;
    }

    typedef GenericGeometry::Pyramid< BaseTopology > Topology;

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

    unsigned int operator() ( const unsigned int order ) const
    {
      return numBaseFunctions_[ order ];
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



  // MonomialBasisHelper
  // -------------------


  template< int mydim, int dim, class F >
  struct MonomialBasisHelper
  {
    typedef MonomialBasisSize< typename GenericGeometry::SimplexTopology< mydim >::type > MySize;
    typedef MonomialBasisSize< typename GenericGeometry::SimplexTopology< dim >::type > Size;

    static void copy ( const unsigned int deriv, F *&wit, F *&rit,
                       const unsigned int numBaseFunctions, const F &z )
    {
      // n(d,k) = size<k>[d];
      MySize &mySize = MySize::instance();
      Size &size = Size::instance();

      const F *const rend = rit + size( deriv )*numBaseFunctions;
      for( ; rit != rend; )
      {
        F *prit = rit;

        *wit = z * *rit;
        ++rit, ++wit;

        for( unsigned d = 1; d <= deriv; ++d )
        {
          #ifndef NDEBUG
          const F *const derivEnd = rit + mySize.sizes_[ d ];
          #endif
          const F *const drend = rit + mySize.sizes_[ d ] - mySize.sizes_[ d-1 ];
          for( ; rit != drend ; ++rit, ++wit )
            *wit = z * *rit;
          for (unsigned int j=1; j<d; ++j)
          {
            const F *const drend = rit + mySize.sizes_[ d-j ] - mySize.sizes_[ d-j-1 ];
            for( ; rit != drend ; ++prit, ++rit, ++wit )
              *wit = F(j) * *prit + z * *rit;
          }
          *wit = F(d) * *prit + z * *rit;
          ++prit, ++rit, ++wit;
          assert(derivEnd == rit);
          rit += size.sizes_[d] - mySize.sizes_[d];
          prit += size.sizes_[d-1] - mySize.sizes_[d-1];
          const F *const emptyWitEnd = wit + size.sizes_[d] - mySize.sizes_[d];
          for ( ; wit != emptyWitEnd; ++wit )
            *wit = Zero<F>();
        }
      }
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

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int block, const unsigned int *const offsets,
                    Field *const values ) const
    {
      *values = Unity< F >();
      F *const end = values + block;
      for( Field *it = values+1 ; it != end; ++it )
        *it = Zero< F >();
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

    typedef MonomialBasisSize< BaseTopology > BaseSize;
    typedef MonomialBasisSize< Topology > Size;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;

    MonomialBasisImpl ()
    {}

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int block, const unsigned int *const offsets,
                    Field *const values ) const
    {
      typedef MonomialBasisHelper< dimDomain, dimD, Field > Helper;
      const BaseSize &size = BaseSize::instance();

      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( deriv, order, x, block, offsets, values );

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *row1 = values + block*offsets[ k-1 ];
        Field *wit = row1 + block*size.sizes_[ k ];
        Helper::copy( deriv, wit, row1, k*size.sizes_[ k ], z );
        Helper::copy( deriv, wit, row0, size( k-1 ), z );
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

    typedef MonomialBasisSize< BaseTopology > BaseSize;
    typedef MonomialBasisSize< Topology > Size;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;

    MonomialBasisImpl ()
    {}

    template< int dimD >
    void evaluateSimplexBase ( const unsigned int deriv, const unsigned int order,
                               const FieldVector< Field, dimD > &x,
                               const unsigned int block, const unsigned int *const offsets,
                               Field *const values,
                               const BaseSize &size ) const
    {
      baseBasis_.evaluate( deriv, order, x, block, offsets, values );
    }

    template< int dimD >
    void evaluatePyramidBase ( const unsigned int deriv, const unsigned int order,
                               const FieldVector< Field, dimD > &x,
                               const unsigned int block, const unsigned int *const offsets,
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
        baseBasis_.evaluate( deriv, order, y, block, offsets, values );

        Field omzk = omz;
        for( unsigned int k = 1; k <= order; ++k )
        {
          Field *it = values + block*offsets[ k-1 ];
          Field *const end = it + block*size.sizes_[ k ];
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
          Field *it = values + block*offsets[ k-1 ];
          Field *const end = it + block*size.sizes_[ k ];
          for( ; it != end; ++it )
            *it = Zero< Field >();
        }
      }
    }

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int block, const unsigned int *const offsets,
                    Field *const values ) const
    {
      typedef MonomialBasisHelper< dimDomain, dimD, Field > Helper;
      const BaseSize &size = BaseSize::instance();

      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplexBase( deriv, order, x, block, offsets, values, size );
      else
        evaluatePyramidBase( deriv, order, x, block, offsets, values, size );

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *row1 = values + block*offsets[ k-1 ];
        Field *wit = row1 + block*size.sizes_[ k ];
        Helper::copy( deriv, wit, row0, size( k-1 ), x[ dimDomain-1 ] );
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
    static const unsigned int dimRange = 1;

    typedef typename Base::Field Field;

    typedef typename Base::DomainVector DomainVector;

    typedef Dune::FieldVector<Field,dimRange> RangeVector;

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

    const unsigned int *sizes () const
    {
      return sizes( order_ );
    }

    const unsigned int size () const
    {
      size_.computeSizes( order_ );
      return size_( order_ );
    }

    const unsigned int derivSize ( const unsigned int deriv ) const
    {
      typedef typename GenericGeometry::SimplexTopology< dimension >::type SimplexTopology;
      MonomialBasisSize< SimplexTopology >::instance().computeSizes( deriv );
      return MonomialBasisSize< SimplexTopology >::instance() ( deriv );
    }

    const unsigned int order () const
    {
      return order_ ;
    }

    void evaluate ( const unsigned int deriv, const DomainVector &x,
                    Field *const values ) const
    {
      Base::evaluate( deriv, order_, x, derivSize( deriv ), sizes( order_ ), values );
    }

    template <unsigned int deriv>
    void evaluate ( const DomainVector &x,
                    Field *const values ) const
    {
      evaluate( deriv, x, values );
    }

    template<unsigned int deriv, class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<deriv>(x,&(values[0]));
    }
    template<unsigned int deriv, DerivativeLayout layout >
    void evaluate ( const DomainVector &x,
                    Derivatives<Field,dimension,1,deriv,layout> *values ) const
    {
      evaluate<deriv>(x,&(values->block()));
    }
    template< unsigned int deriv >
    void evaluate ( const DomainVector &x,
                    FieldVector<Field,Derivatives<Field,dimension,1,deriv,value>::size> *values ) const
    {
      evaluate(0,x,&(values[0][0]));
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
    template <class Vector>
    void integral ( Vector &values ) const
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
    static const unsigned int dimRange = 1;

    typedef FieldVector<Field,dimension> DomainVector;
    typedef FieldVector<Field,dimRange> RangeVector;

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
    template < unsigned int deriv >
    void evaluate ( const DomainVector &x,
                    Field *const values ) const
    {
      evaluate( deriv, x, values );
    }
    template < unsigned int deriv, int size >
    void evaluate ( const DomainVector &x,
                    Dune::FieldVector<Field,size> *const values ) const
    {
      evaluate( deriv, x, &(values[0][0]) );
    }
    template<unsigned int deriv, DerivativeLayout layout >
    void evaluate ( const DomainVector &x,
                    Derivatives<Field,dimension,1,deriv,layout> *values ) const
    {
      evaluate<deriv>(x,&(values->block()));
    }
    template <unsigned int deriv, class Vector>
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<deriv>( x, &(values[ 0 ]) );
    }
    template< class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<0>(x,values);
    }

    virtual void integral ( Field *const values ) const = 0;
    template <class Vector>
    void integral ( Vector &values ) const
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



  // MonomialBasisCreator
  // --------------------

  template< int dim, class F >
  struct MonomialBasisCreator
  {
    static const unsigned int dimension = dim;

    typedef F StorageField;

    typedef unsigned int Key;
    typedef VirtualMonomialBasis< dimension, StorageField > Basis;

    template< class Topology >
    static const Basis &basis ( const Key &order )
    {
      return *(new VirtualMonomialBasisImpl< Topology, StorageField >( order ));
    }

    static void release ( const Basis &basis )
    {
      delete &basis;
    }

    template< class Topology >
    static void basis(unsigned int order,Basis* &basis)
    {
      basis = new VirtualMonomialBasisImpl<Topology,StorageField>(order);
    }
  };



  // MonomialBasisProvider
  // ---------------------

  template< int dim, class SF >
  struct MonomialBasisProvider
    : public BasisProvider< MonomialBasisCreator< dim, SF > >
  {};

}

#endif
