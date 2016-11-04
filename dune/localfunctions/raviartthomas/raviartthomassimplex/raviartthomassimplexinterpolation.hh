// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH

#include <fstream>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/utility/interpolationhelper.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>
#include <dune/localfunctions/orthonormal/orthonormalbasis.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactory;

  template < unsigned int dim, class Field >
  struct RaviartThomasL2InterpolationFactory;



  // LocalCoefficientsContainer
  // --------------------------

  class LocalCoefficientsContainer
  {
    typedef LocalCoefficientsContainer This;

  public:
    template <class Setter>
    LocalCoefficientsContainer ( const Setter &setter )
    {
      setter.setLocalKeys(localKey_);
    }

    const LocalKey &localKey ( const unsigned int i ) const
    {
      assert( i < size() );
      return localKey_[ i ];
    }

    unsigned int size () const
    {
      return localKey_.size();
    }

  private:
    std::vector< LocalKey > localKey_;
  };



  // RaviartThomasCoefficientsFactoryTraits
  // --------------------------------------

  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactoryTraits
  {
    static const unsigned int dimension = dim;
    typedef const LocalCoefficientsContainer Object;
    typedef unsigned int Key;
    typedef RaviartThomasCoefficientsFactory<dim> Factory;
  };



  // RaviartThomasCoefficientsFactory
  // --------------------------------

  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactory
    : public TopologyFactory< RaviartThomasCoefficientsFactoryTraits< dim > >
  {
    typedef RaviartThomasCoefficientsFactoryTraits< dim > Traits;

    template< class Topology >
    static typename Traits::Object *createObject( const typename Traits::Key &key )
    {
      typedef RaviartThomasL2InterpolationFactory< dim, double > InterpolationFactory;
      if( !supports< Topology >( key ) )
        return nullptr;
      typename InterpolationFactory::Object *interpolation = InterpolationFactory::template create< Topology >( key );
      typename Traits::Object *localKeys = new typename Traits::Object( *interpolation );
      InterpolationFactory::release( interpolation );
      return localKeys;
    }

    template< class Topology >
    static bool supports ( const typename Traits::Key &key )
    {
      return Impl::IsSimplex< Topology >::value;
    }
  };



  // RTL2InterpolationBuilder
  // ------------------------

  // L2 Interpolation requires:
  // - for element
  //   - test basis
  // - for each face (dynamic)
  //   - test basis
  //   - normal
  template< unsigned int dim, class Field >
  struct RTL2InterpolationBuilder
  {
    static const unsigned int dimension = dim;
    typedef OrthonormalBasisFactory< dimension, Field > TestBasisFactory;
    typedef typename TestBasisFactory::Object TestBasis;
    typedef OrthonormalBasisFactory< dimension-1, Field > TestFaceBasisFactory;
    typedef typename TestFaceBasisFactory::Object TestFaceBasis;
    typedef FieldVector< Field, dimension > Normal;

    RTL2InterpolationBuilder () = default;

    RTL2InterpolationBuilder ( const RTL2InterpolationBuilder & ) = delete;
    RTL2InterpolationBuilder ( RTL2InterpolationBuilder && ) = delete;

    ~RTL2InterpolationBuilder ()
    {
      TestBasisFactory::release( testBasis_ );
      for( FaceStructure &f : faceStructure_ )
        TestFaceBasisFactory::release( f.basis_ );
    }

    unsigned int topologyId () const { return topologyId_; }

    GeometryType type () const { return GeometryType( topologyId(), dimension ); }

    unsigned int order () const { return order_; }

    unsigned int faceSize () const { return faceSize_; }

    TestBasis *testBasis () const { return testBasis_; }
    TestFaceBasis *testFaceBasis ( unsigned int f ) const { assert( f < faceSize() ); return faceStructure_[ f ].basis_; }

    const Normal &normal ( unsigned int f ) const { assert( f < faceSize() ); return *(faceStructure_[ f ].normal_); }

    template< class Topology >
    void build ( unsigned int order )
    {
      order_ = order;
      topologyId_ = Topology::id;

      testBasis_ = (order > 0 ? TestBasisFactory::template create< Topology >( order-1 ) : nullptr);

      const auto &refElement = ReferenceElements< Field, dimension >::general( type() );
      faceSize_ = refElement.size( 1 );
      faceStructure_.reserve( faceSize_ );
      for( unsigned int face = 0; face < faceSize_; ++face )
      {
        TestFaceBasis *faceBasis = Impl::IfTopology< CreateFaceBasis, dimension-1 >::apply( refElement.type( face, 1 ).id(), order );
        faceStructure_.emplace_back( faceBasis, refElement.integrationOuterNormal( face ) );
      }
      assert( faceStructure_.size() == faceSize_ );
    }

  private:
    struct FaceStructure
    {
      FaceStructure( TestFaceBasis *tfb, const Normal &n )
        : basis_( tfb ), normal_( &n )
      {}

      TestFaceBasis *basis_;
      const Dune::FieldVector< Field, dimension > *normal_;
    };

    template< class FaceTopology >
    struct CreateFaceBasis
    {
      static TestFaceBasis *apply ( unsigned int order ) { return TestFaceBasisFactory::template create< FaceTopology >( order ); }
    };

    std::vector< FaceStructure > faceStructure_;
    TestBasis *testBasis_ = nullptr;
    unsigned int topologyId_, order_, faceSize_;
  };



  // RaviartThomasL2Interpolation
  // ----------------------------

  /**
   * \class RaviartThomasL2Interpolation
   * \brief An L2-based interpolation for Raviart Thomas
   *
   **/
  template< unsigned int dimension, class F>
  class RaviartThomasL2Interpolation
    : public InterpolationHelper< F ,dimension >
  {
    typedef RaviartThomasL2Interpolation< dimension, F > This;
    typedef InterpolationHelper<F,dimension> Base;

  public:
    typedef F Field;
    typedef RTL2InterpolationBuilder<dimension,Field> Builder;
    RaviartThomasL2Interpolation()
      : order_(0),
        size_(0)
    {}

    template< class Function, class Fy >
    void interpolate ( const Function &function, std::vector< Fy > &coefficients ) const
    {
      coefficients.resize(size());
      typename Base::template Helper<Function,std::vector<Fy>,true> func( function,coefficients );
      interpolate(func);
    }
    template< class Basis, class Matrix >
    void interpolate ( const Basis &basis, Matrix &matrix ) const
    {
      matrix.resize( size(), basis.size() );
      typename Base::template Helper<Basis,Matrix,false> func( basis,matrix );
      interpolate(func);
    }

    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return size_;
    }
    template <class Topology>
    void build( unsigned int order )
    {
      size_ = 0;
      order_ = order;
      builder_.template build<Topology>(order_);
      if (builder_.testBasis())
        size_ += dimension*builder_.testBasis()->size();
      for ( unsigned int f=0; f<builder_.faceSize(); ++f )
        if (builder_.testFaceBasis(f))
          size_ += builder_.testFaceBasis(f)->size();
    }

    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      unsigned int row = 0;
      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        if (builder_.faceSize())
          for (unsigned int i=0; i<builder_.testFaceBasis(f)->size(); ++i,++row)
            keys[row] = LocalKey(f,1,i);
      }
      if (builder_.testBasis())
        for (unsigned int i=0; i<builder_.testBasis()->size()*dimension; ++i,++row)
          keys[row] = LocalKey(0,0,i);
      assert( row == size() );
    }

  protected:
    template< class Func, class Container, bool type >
    void interpolate ( typename Base::template Helper<Func,Container,type> &func ) const
    {
      const Dune::GeometryType geoType( builder_.topologyId(), dimension );

      std::vector< Field > testBasisVal;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // boundary dofs:
      typedef Dune::QuadratureRule<Field, dimension-1> FaceQuadrature;
      typedef Dune::QuadratureRules<Field, dimension-1> FaceQuadratureRules;

      typedef Dune::ReferenceElements< Field, dimension > RefElements;
      typedef Dune::ReferenceElement< Field, dimension > RefElement;
      typedef typename RefElement::template Codim< 1 >::Geometry Geometry;

      const RefElement &refElement = RefElements::general( geoType );

      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        if (!builder_.testFaceBasis(f))
          continue;
        testBasisVal.resize(builder_.testFaceBasis(f)->size());

        const Geometry &geometry = refElement.template geometry< 1 >( f );
        const Dune::GeometryType subGeoType( geometry.type().id(), dimension-1 );
        const FaceQuadrature &faceQuad = FaceQuadratureRules::rule( subGeoType, 2*order_+2 );

        const unsigned int quadratureSize = faceQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          if (dimension>1)
            builder_.testFaceBasis(f)->template evaluate<0>(faceQuad[qi].position(),testBasisVal);
          else
            testBasisVal[0] = 1.;
          fillBnd( row, testBasisVal,
                   func.evaluate( geometry.global( faceQuad[qi].position() ) ),
                   builder_.normal(f), faceQuad[qi].weight(),
                   func);
        }

        row += builder_.testFaceBasis(f)->size();
      }
      // element dofs
      if (builder_.testBasis())
      {
        testBasisVal.resize(builder_.testBasis()->size());

        typedef Dune::QuadratureRule<Field, dimension> Quadrature;
        typedef Dune::QuadratureRules<Field, dimension> QuadratureRules;
        const Quadrature &elemQuad = QuadratureRules::rule( geoType, 2*order_+1 );

        const unsigned int quadratureSize = elemQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          builder_.testBasis()->template evaluate<0>(elemQuad[qi].position(),testBasisVal);
          fillInterior( row, testBasisVal,
                        func.evaluate(elemQuad[qi].position()),
                        elemQuad[qi].weight(),
                        func );
        }

        row += builder_.testBasis()->size()*dimension;
      }
      assert(row==size());
    }

  private:
    /** /brief evaluate boundary functionals **/
    template <class MVal, class RTVal,class Matrix>
    void fillBnd (unsigned int startRow,
                  const MVal &mVal,
                  const RTVal &rtVal,
                  const FieldVector<Field,dimension> &normal,
                  const Field &weight,
                  Matrix &matrix) const
    {
      const unsigned int endRow = startRow+mVal.size();
      typename RTVal::const_iterator rtiter = rtVal.begin();
      for ( unsigned int col = 0; col < rtVal.size() ; ++rtiter,++col)
      {
        Field cFactor = (*rtiter)*normal;
        typename MVal::const_iterator miter = mVal.begin();
        for (unsigned int row = startRow;
             row!=endRow; ++miter, ++row )
        {
          matrix.add(row,col, (weight*cFactor)*(*miter) );
        }
        assert( miter == mVal.end() );
      }
    }
    template <class MVal, class RTVal,class Matrix>
    void fillInterior (unsigned int startRow,
                       const MVal &mVal,
                       const RTVal &rtVal,
                       Field weight,
                       Matrix &matrix) const
    {
      const unsigned int endRow = startRow+mVal.size()*dimension;
      typename RTVal::const_iterator rtiter = rtVal.begin();
      for ( unsigned int col = 0; col < rtVal.size() ; ++rtiter,++col)
      {
        typename MVal::const_iterator miter = mVal.begin();
        for (unsigned int row = startRow;
             row!=endRow; ++miter,row+=dimension )
        {
          for (unsigned int i=0; i<dimension; ++i)
          {
            matrix.add(row+i,col, (weight*(*miter))*(*rtiter)[i] );
          }
        }
        assert( miter == mVal.end() );
      }
    }

    Builder builder_;
    unsigned int order_;
    unsigned int size_;
  };

  template < unsigned int dim, class F >
  struct RaviartThomasL2InterpolationFactoryTraits
  {
    static const unsigned int dimension = dim;
    typedef unsigned int Key;
    typedef const RaviartThomasL2Interpolation<dim,F> Object;
    typedef RaviartThomasL2InterpolationFactory<dim,F> Factory;
  };
  template < unsigned int dim, class Field >
  struct RaviartThomasL2InterpolationFactory :
    public TopologyFactory< RaviartThomasL2InterpolationFactoryTraits<dim,Field> >
  {
    typedef RaviartThomasL2InterpolationFactoryTraits<dim,Field> Traits;
    typedef RTL2InterpolationBuilder<dim,Field> Builder;
    typedef typename Traits::Object Object;
    typedef typename std::remove_const<Object>::type NonConstObject;
    template <class Topology>
    static typename Traits::Object *createObject( const typename Traits::Key &key )
    {
      if ( !supports<Topology>(key) )
        return 0;
      NonConstObject *interpol = new NonConstObject();
      interpol->template build<Topology>(key);
      return interpol;
    }
    template< class Topology >
    static bool supports ( const typename Traits::Key &key )
    {
      return Impl::IsSimplex<Topology>::value;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH
