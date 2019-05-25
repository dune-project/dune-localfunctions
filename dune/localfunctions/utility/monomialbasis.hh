// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/multiindex.hh>
#include <dune/localfunctions/utility/tensor.hh>

namespace Dune
{
  /************************************************
  * Classes for evaluating ''Monomials'' on any order
  * for all reference element type.
  * For a simplex topology these are the normal
  * monomials for cube topologies the bimonomials.
  * The construction follows the construction of the
  * generic geometries using tensor products for
  * prism generation and duffy transform for pyramid
  * construction.
  * A derivative argument can be applied, in which case
  * all derivatives up to the desired order are
  * evaluated. Note that for higher order derivatives
  * only the ''lower'' part of the symmetric tensor
  * is evaluated, e.g., passing derivative equal to 2
  * to the class will provide the vector
  *    (d/dxdx p, d/dxydx p, d/dydy p,
  *     d/dx p, d/dy p, p)
  * Important:
  * So far the computation of the derivatives has not
  * been fully implemented for general pyramid
  * construction, i.e., in the case where a pyramid is
  * build over a non simplex base geometry.
  *
  * Central classes:
  * 1) template< class Topology, class F >
  *    class MonomialBasisImpl;
  *    Implementation of the monomial evaluation for
  *    a given topology and field type.
  *    The method evaluate fills a F* vector
  * 2) template< class Topology, class F >
  *    class MonomialBasis
  *    The base class for the static monomial evaluation
  *    providing addiional evaluate methods including
  *    one taking std::vector<F>.
  * 3) template< int dim, class F >
  *    class VirtualMonomialBasis
  *    Virtualization of the MonomialBasis.
  * 4) template< int dim, class F >
  *    struct MonomialBasisFactory;
  *    A factory class for the VirtualMonomialBasis
  * 5) template< int dim, class F >
  *    struct MonomialBasisProvider
  *    A singleton container for the virtual monomial
  *    basis
  ************************************************/

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology >
  class MonomialBasisSize;

  template< class Topology, class F >
  class MonomialBasis;



  // MonomialBasisSize
  // -----------------

  template< class TopologyType >
  class MonomialBasisSize
  {
    typedef MonomialBasisSize< TopologyType > This;

  public:
    static This &instance ()
    {
      static This _instance;
      return _instance;
    }

    unsigned int maxOrder_;

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

    void computeSizes ( unsigned int order )
    {
      if (order <= maxOrder_)
        return;

      maxOrder_ = order;

      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      constexpr auto dim = TopologyType::dimension;

      sizes_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
        sizes_[ k ] = 0;

      std::fill(numBaseFunctions_, numBaseFunctions_+order+1, 1);

      for( int codim=dim-1; codim>=0; codim--)
      {
        if (Impl::isPrism(TopologyType::id,dim,codim))
        {
          for( unsigned int k = 1; k <= order; ++k )
          {
            sizes_[ k ]            = numBaseFunctions_[ k ] + k*sizes_[ k ];
            numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
          }
        }
        else
        {
          for( unsigned int k = 1; k <= order; ++k )
          {
            sizes_[ k ]            = numBaseFunctions_[ k ];
            numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
          }
        }
      }
    }
  };



  // MonomialBasisHelper
  // -------------------


  template< int mydim, int dim, class F >
  struct MonomialBasisHelper
  {
    typedef MonomialBasisSize< typename Impl::SimplexTopology< mydim >::type > MySize;
    typedef MonomialBasisSize< typename Impl::SimplexTopology< dim >::type > Size;

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

          {
            const F *const drend = rit + mySize.sizes_[ d ] - mySize.sizes_[ d-1 ];
            for( ; rit != drend ; ++rit, ++wit )
              *wit = z * *rit;
          }

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
  class MonomialBasisImpl< Impl::Point, F >
  {
    typedef MonomialBasisImpl< Impl::Point, F > This;

  public:
    typedef Impl::Point Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< Impl::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< Impl::Pyramid< Topology >, Field >;

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

    void integrate ( const unsigned int order,
                     const unsigned int *const offsets,
                     Field *const values ) const
    {
      values[ 0 ] = Unity< Field >();
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasisImpl< Impl::Prism< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< Impl::Prism< BaseTopology >, F > This;

  public:
    typedef Impl::Prism< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< Impl::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< Impl::Pyramid< Topology >, Field >;

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
      const_cast<BaseSize&>(size).computeSizes(order);

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

    void integrate ( const unsigned int order,
                     const unsigned int *const offsets,
                     Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();
      const Size &mySize = Size::instance();
      // fill first column
      baseBasis_.integrate( order, offsets, values );
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
  class MonomialBasisImpl< Impl::Pyramid< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< Impl::Pyramid< BaseTopology >, F > This;

  public:
    typedef Impl::Pyramid< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< Impl::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< Impl::Pyramid< Topology >, Field >;

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
      const_cast<BaseSize&>(size).computeSizes(order);

      if( Impl::IsSimplex< Topology >::value )
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

    void integrate ( const unsigned int order,
                     const unsigned int *const offsets,
                     Field *const values ) const
    {
      const BaseSize &size = BaseSize::instance();

      // fill first column
      baseBasis_.integrate( order, offsets, values );

      const unsigned int *const baseSizes = size.sizes_;

      {
        Field *const col0End = values + baseSizes[ 0 ];
        for( Field *it = values; it != col0End; ++it )
          *it *= Field( 1 ) /  Field( int(dimDomain) );
      }

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
    {
      assert(order<=1024); // avoid wrapping of unsigned int (0-1) order=1024 is quite hight...)
    }

    const unsigned int *sizes ( unsigned int order ) const
    {
      size_.computeSizes( order );
      return size_.numBaseFunctions_;
    }

    const unsigned int *sizes () const
    {
      return sizes( order_ );
    }

    unsigned int size () const
    {
      size_.computeSizes( order_ );
      return size_( order_ );
    }

    unsigned int derivSize ( const unsigned int deriv ) const
    {
      typedef typename Impl::SimplexTopology< dimension >::type SimplexTopology;
      MonomialBasisSize< SimplexTopology >::instance().computeSizes( deriv );
      return MonomialBasisSize< SimplexTopology >::instance() ( deriv );
    }

    unsigned int order () const
    {
      return order_ ;
    }

    unsigned int topologyId ( ) const
    {
      return Topology::id;
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
    template<unsigned int deriv, DerivativeLayoutNS::DerivativeLayout layout >
    void evaluate ( const DomainVector &x,
                    Derivatives<Field,dimension,1,deriv,layout> *values ) const
    {
      evaluate<deriv>(x,&(values->block()));
    }
    template< unsigned int deriv >
    void evaluate ( const DomainVector &x,
                    FieldVector<Field,Derivatives<Field,dimension,1,deriv,DerivativeLayoutNS::value>::size> *values ) const
    {
      evaluate(0,x,&(values[0][0]));
    }

    template<class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<0>(x,&(values[0]));
    }

    template< class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<0>( bx, values );
    }

    void integrate ( Field *const values ) const
    {
      Base::integrate( order_, sizes( order_ ), values );
    }
    template <class Vector>
    void integrate ( Vector &values ) const
    {
      integrate( &(values[ 0 ]) );
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
    : public MonomialBasis< typename Impl::SimplexTopology< dim >::type, F >
  {
    typedef StandardMonomialBasis< dim, F > This;
    typedef MonomialBasis< typename Impl::SimplexTopology< dim >::type, F > Base;

  public:
    typedef typename Impl::SimplexTopology< dim >::type Topology;
    static const int dimension = dim;

    StandardMonomialBasis ( unsigned int order )
      : Base( order )
    {}
  };



  // StandardBiMonomialBasis
  // -----------------------

  template< int dim, class F >
  class StandardBiMonomialBasis
    : public MonomialBasis< typename Impl::CubeTopology< dim >::type, F >
  {
    typedef StandardBiMonomialBasis< dim, F > This;
    typedef MonomialBasis< typename Impl::CubeTopology< dim >::type, F > Base;

  public:
    typedef typename Impl::CubeTopology< dim >::type Topology;
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
  class VirtualMonomialBasis
  {
    typedef VirtualMonomialBasis< dim, F > This;

  public:
    typedef F Field;
    typedef F StorageField;
    static const int dimension = dim;
    static const unsigned int dimRange = 1;

    typedef FieldVector<Field,dimension> DomainVector;
    typedef FieldVector<Field,dimRange> RangeVector;

    explicit VirtualMonomialBasis(unsigned int topologyId,
                                  unsigned int order)
      : order_(order), topologyId_(topologyId) {}

    virtual ~VirtualMonomialBasis() {}

    virtual const unsigned int *sizes ( ) const = 0;

    unsigned int size ( ) const
    {
      return sizes( )[ order_ ];
    }

    unsigned int order () const
    {
      return order_;
    }

    unsigned int topologyId ( ) const
    {
      return topologyId_;
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
    template<unsigned int deriv, DerivativeLayoutNS::DerivativeLayout layout >
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
    template< class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<0>( bx, values );
    }
    template< unsigned int deriv, class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<deriv>( bx, values );
    }

    virtual void integrate ( Field *const values ) const = 0;
    template <class Vector>
    void integrate ( Vector &values ) const
    {
      integrate( &(values[ 0 ]) );
    }
  protected:
    unsigned int order_;
    unsigned int topologyId_;
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
      : Base(Topology::id,order), basis_(order)
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

    void integrate ( Field *const values ) const
    {
      basis_.integrate(values);
    }

  private:
    MonomialBasis<Topology,Field> basis_;
    using Base::order_;
  };

  // MonomialBasisFactory
  // --------------------

  template< int dim, class F >
  struct MonomialBasisFactory
  {
    static const unsigned int dimension = dim;
    typedef F StorageField;

    typedef unsigned int Key;
    typedef const VirtualMonomialBasis< dimension, F > Object;

    template < int dd, class FF >
    struct EvaluationBasisFactory
    {
      typedef MonomialBasisFactory<dd,FF> Type;
    };

    template< class Topology >
    static Object* create ( const Key &order )
    {
      return new VirtualMonomialBasisImpl< Topology, StorageField >( order );
    }
    static void release( Object *object ) { delete object; }
  };



  // MonomialBasisProvider
  // ---------------------

  template< int dim, class SF >
  struct MonomialBasisProvider
    : public TopologySingletonFactory< MonomialBasisFactory< dim, SF > >
  {
    static const unsigned int dimension = dim;
    typedef SF StorageField;
    template < int dd, class FF >
    struct EvaluationBasisFactory
    {
      typedef MonomialBasisProvider<dd,FF> Type;
    };
  };

}

#endif
