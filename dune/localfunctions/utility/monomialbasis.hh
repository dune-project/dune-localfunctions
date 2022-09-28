// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/topologyfactory.hh>

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
  * 1) template< GeometryType::Id geometryId, class F >
  *    class MonomialBasisImpl;
  *    Implementation of the monomial evaluation for
  *    a given topology and field type.
  *    The method evaluate fills a F* vector
  * 2) template< GeometryType::Id geometryId, class F >
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

  template< GeometryType::Id geometryId >
  class MonomialBasisSize;

  template< GeometryType::Id geometryId, class F >
  class MonomialBasis;



  // MonomialBasisSize
  // -----------------

  template< GeometryType::Id geometryId >
  class MonomialBasisSize
  {
    typedef MonomialBasisSize< geometryId > This;

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

      constexpr GeometryType geometry = geometryId;
      constexpr auto dim = geometry.dim();

      sizes_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
        sizes_[ k ] = 0;

      std::fill(numBaseFunctions_, numBaseFunctions_+order+1, 1);

      for( int codim=dim-1; codim>=0; codim--)
      {
        if (Impl::isPrism(geometry.id(),dim,codim))
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
    typedef MonomialBasisSize< GeometryTypes::simplex(mydim).toId() > MySize;
    typedef MonomialBasisSize< GeometryTypes::simplex(dim).toId() > Size;

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

  template< GeometryType::Id geometryId, class F>
  class MonomialBasisImpl
  {
  public:
    typedef F Field;

    static constexpr GeometryType geometry = geometryId;

    static const unsigned int dimDomain = geometry.dim();

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< geometryId, Field >;

    MonomialBasisImpl ()
    {}

    template< int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int block, const unsigned int *const offsets,
                    Field *const values ) const
    {
      //start with vertex
      *values = Unity< F >();
      F *const end = values + block;
      for( Field *it = values+1 ; it != end; ++it )
        *it = Zero< F >();

      constexpr GeometryType gt = GeometryTypes::vertex;

      if constexpr ( geometry == gt)
        return;
      else
        evaluate<gt,dimD>(deriv, order, x, block, offsets, values );
    }

    template<GeometryType::Id baseGeometryId, int dimD >
    void evaluate ( const unsigned int deriv, const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int block, const unsigned int *const offsets,
                    Field *const values ) const
    {

      static constexpr GeometryType baseGeometry = baseGeometryId;

      auto constexpr isPrismatic = geometry.isPrismatic(baseGeometry.dim());

      // compute
      typedef MonomialBasisHelper< baseGeometry.dim() + 1, dimD, Field > Helper;
      typedef MonomialBasisSize<baseGeometryId> BaseSize;

      const BaseSize &size = BaseSize::instance();
      const_cast<BaseSize&>(size).computeSizes(order);

      const Field &z = x[ baseGeometry.dim() ];

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *row1 = values + block*offsets[ k-1 ];
        Field *wit = row1 + block*size.sizes_[ k ];
        if constexpr ( isPrismatic )
          Helper::copy( deriv, wit, row1, k*size.sizes_[ k ], z );
        Helper::copy( deriv, wit, row0, size( k-1 ), z );
        row0 = row1;
      }

      // stop if desired dimension is reached
      if constexpr( baseGeometry.dim() == dimDomain-1)
        return;
      else
      {
        constexpr GeometryType nextGeometry = isPrismatic ? GeometryTypes::prismaticExtension(baseGeometry)
                                                          : GeometryTypes::conicalExtension(baseGeometry);

        evaluate<nextGeometry.toId(),dimD>(deriv, order, x, block, offsets, values );
      }
    }

    void integrate ( const unsigned int order,
                     const unsigned int *const offsets,
                     Field *const values ) const
    {
      //start with vertex
      values[ 0 ] = Unity< Field >();
      static constexpr GeometryType gt = GeometryTypes::vertex;

      if constexpr ( geometry == gt)
        return;
      else
        integrate<gt>(order, offsets, values);
    }

    template<GeometryType::Id baseGeometryId>
    void integrate ( const unsigned int order,
                     const unsigned int *const offsets,
                     Field *const values) const
    {
      static constexpr GeometryType baseGeometry = baseGeometryId;

      auto constexpr isPrismatic = geometry.isPrismatic(baseGeometry.dim());

      // decide which kind of integration should be performed
      if constexpr ( isPrismatic )
        integratePrismatic<baseGeometry>(order, offsets, values);
      else
        integrateConical<baseGeometry>(order, offsets, values);

      // stop if the desired dimension is reached
      if constexpr( baseGeometry.dim() == dimDomain-1)
        return;
      else
      {
        static constexpr GeometryType nextGeometry = (isPrismatic ? GeometryTypes::prismaticExtension(baseGeometry)
                                                                  : GeometryTypes::conicalExtension(baseGeometry));

        integrate<nextGeometry.toId()>(order, offsets, values);
      }

    }

    template<GeometryType::Id baseGeometryId>
    void integratePrismatic ( const unsigned int order,
                              const unsigned int *const offsets,
                              Field *const values ) const
    {
      typedef MonomialBasisSize<baseGeometryId> BaseSize;
      static const BaseSize &size = BaseSize::instance();
      const unsigned int *const baseSizes = size.sizes_;

      static constexpr GeometryType baseGeometry = baseGeometryId;
      static constexpr GeometryType nextGeometry = GeometryTypes::prismaticExtension(baseGeometry);

      typedef MonomialBasisSize<nextGeometry.toId()> Size;
      static const Size &mySize = Size::instance();

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


    template<GeometryType::Id baseGeometryId>
    void integrateConical ( const unsigned int order,
                            const unsigned int *const offsets,
                            Field *const values) const
    {
      typedef MonomialBasisSize<baseGeometryId> BaseSize;
      static const BaseSize &size = BaseSize::instance();
      const unsigned int *const baseSizes = size.sizes_;

      static constexpr GeometryType baseGeometry = baseGeometryId;

      {
        Field *const col0End = values + baseSizes[ 0 ];
        for( Field *it = values; it != col0End; ++it )
          *it *= Field( 1 ) /  Field( int(baseGeometry.dim()+1) );
      }

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        const Field factor = (Field( 1 ) / Field( k + baseGeometry.dim()+1));

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

  template< GeometryType::Id geometryId, class F >
  class MonomialBasis
    : public MonomialBasisImpl< geometryId, F >
  {
    static constexpr GeometryType geometry = geometryId;
    typedef MonomialBasis< geometryId, F > This;
    typedef MonomialBasisImpl< geometryId, F > Base;

  public:
    static const unsigned int dimension = Base::dimDomain;
    static const unsigned int dimRange = 1;

    typedef typename Base::Field Field;

    typedef typename Base::DomainVector DomainVector;

    typedef Dune::FieldVector<Field,dimRange> RangeVector;

    typedef MonomialBasisSize<geometryId> Size;

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
      MonomialBasisSize< GeometryTypes::simplex(dimension).toId() >::instance().computeSizes( deriv );
      return MonomialBasisSize< GeometryTypes::simplex(dimension).toId() >::instance() ( deriv );
    }

    unsigned int order () const
    {
      return order_ ;
    }

    unsigned int topologyId ( ) const
    {
      return geometry.id();
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
    : public MonomialBasis< GeometryTypes::simplex(dim).toId() , F >
  {
    typedef StandardMonomialBasis< dim, F > This;
    typedef MonomialBasis< GeometryTypes::simplex(dim).toId(), F > Base;

  public:
    static constexpr GeometryType geometry = GeometryTypes::simplex(dim);
    static const int dimension = dim;

    StandardMonomialBasis ( unsigned int order )
      : Base( order )
    {}
  };



  // StandardBiMonomialBasis
  // -----------------------

  template< int dim, class F >
  class StandardBiMonomialBasis
    : public MonomialBasis< GeometryTypes::cube(dim).toId() , F >
  {
    typedef StandardBiMonomialBasis< dim, F > This;
    typedef MonomialBasis< GeometryTypes::cube(dim).toId() , F > Base;

  public:
    static constexpr GeometryType geometry = GeometryTypes::cube(dim);
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

    explicit VirtualMonomialBasis(const GeometryType& gt,
                                  unsigned int order)
      : order_(order), geometry_(gt) {}

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

    GeometryType type() const
    {
      return geometry_;
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
    GeometryType geometry_;
  };

  template< GeometryType::Id geometryId, class F >
  class VirtualMonomialBasisImpl
    : public VirtualMonomialBasis< static_cast<GeometryType>(geometryId).dim(), F >
  {
    static constexpr GeometryType geometry = geometryId;
    typedef VirtualMonomialBasis< geometry.dim(), F > Base;
    typedef VirtualMonomialBasisImpl< geometryId, F > This;

  public:
    typedef typename Base::Field Field;
    typedef typename Base::DomainVector DomainVector;

    VirtualMonomialBasisImpl(unsigned int order)
      : Base(geometry,order), basis_(order)
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
    MonomialBasis<geometryId,Field> basis_;
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

    template< GeometryType::Id geometryId >
    static Object* create ( const Key &order )
    {
      return new VirtualMonomialBasisImpl< geometryId, StorageField >( order );
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
