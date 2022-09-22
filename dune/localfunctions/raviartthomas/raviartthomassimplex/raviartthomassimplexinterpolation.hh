// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH

#include <fstream>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/utility/interpolationhelper.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>
#include <dune/localfunctions/orthonormal/orthonormalbasis.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

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

    std::size_t size () const
    {
      return localKey_.size();
    }

  private:
    std::vector< LocalKey > localKey_;
  };



  // RaviartThomasCoefficientsFactory
  // --------------------------------

  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactory
  {
    typedef std::size_t Key;
    typedef const LocalCoefficientsContainer Object;

    template< GeometryType::Id geometryId >
    static Object *create( const Key &key )
    {
      typedef RaviartThomasL2InterpolationFactory< dim, double > InterpolationFactory;
      if( !supports< geometryId >( key ) )
        return nullptr;
      typename InterpolationFactory::Object *interpolation = InterpolationFactory::template create< geometryId >( key );
      Object *localKeys = new Object( *interpolation );
      InterpolationFactory::release( interpolation );
      return localKeys;
    }

    template< GeometryType::Id geometryId >
    static bool supports ( const Key &key )
    {
      return GeometryType(geometryId).isSimplex();
    }
    static void release( Object *object ) { delete object; }
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

    // for the dofs associated to the element
    typedef OrthonormalBasisFactory< dimension, Field > TestBasisFactory;
    typedef typename TestBasisFactory::Object TestBasis;

    // for the dofs associated to the faces
    typedef OrthonormalBasisFactory< dimension-1, Field > TestFaceBasisFactory;
    typedef typename TestFaceBasisFactory::Object TestFaceBasis;

    // the normals of the faces
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

    GeometryType type () const { return geometry_; }

    std::size_t order () const { return order_; }

    // number of faces
    unsigned int faceSize () const { return faceSize_; }

    // basis associated to the element
    TestBasis *testBasis () const { return testBasis_; }

    // basis associated to face f
    TestFaceBasis *testFaceBasis ( unsigned int f ) const { assert( f < faceSize() ); return faceStructure_[ f ].basis_; }

    // normal of face f
    const Normal &normal ( unsigned int f ) const { assert( f < faceSize() ); return *(faceStructure_[ f ].normal_); }

    template< GeometryType::Id geometryId >
    void build ( std::size_t order )
    {
      constexpr GeometryType geometry = geometryId;
      geometry_ = geometry;
      order_ = order;

      testBasis_ = (order > 0 ? TestBasisFactory::template create< geometry >( order-1 ) : nullptr);

      const auto &refElement = ReferenceElements< Field, dimension >::general( type() );
      faceSize_ = refElement.size( 1 );
      faceStructure_.reserve( faceSize_ );
      for( unsigned int face = 0; face < faceSize_; ++face )
      {
        /* For simplices or cubes of arbitrary dimension you could just use
         *
         * ```
         * GeometryType faceGeometry = Impl::getBase(geometry_);
         * TestFaceBasis *faceBasis = TestFaceBasisFactory::template create< faceGeometry >( order );
         * ```
         *
         * For i.e. Prisms and Pyramids in 3d this does not work because they contain squares and triangles as faces.
         * And depending on the dynamic face index a different face geometry is needed.
         *
         */
        TestFaceBasis *faceBasis = Impl::toGeometryTypeIdConstant<dimension-1>(refElement.type( face, 1 ), [&](auto faceGeometryTypeId) {
            return TestFaceBasisFactory::template create< decltype(faceGeometryTypeId)::value >( order );
            });
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

    std::vector< FaceStructure > faceStructure_;
    TestBasis *testBasis_ = nullptr;
    GeometryType geometry_;
    unsigned int faceSize_;
    std::size_t order_;
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

    template< class Function, class Vector >
    auto interpolate ( const Function &function, Vector &coefficients ) const
    -> std::enable_if_t< std::is_same< decltype(std::declval<Vector>().resize(1) ),void >::value,void>
    {
      coefficients.resize(size());
      typename Base::template Helper<Function,Vector,true> func( function,coefficients );
      interpolate(func);
    }

    template< class Basis, class Matrix >
    auto interpolate ( const Basis &basis, Matrix &matrix ) const
    -> std::enable_if_t< std::is_same<
           decltype(std::declval<Matrix>().rowPtr(0)), typename Matrix::Field* >::value,void>
    {
      matrix.resize( size(), basis.size() );
      typename Base::template Helper<Basis,Matrix,false> func( basis,matrix );
      interpolate(func);
    }

    std::size_t order() const
    {
      return order_;
    }
    std::size_t size() const
    {
      return size_;
    }
    template <GeometryType::Id geometryId>
    void build( std::size_t order )
    {
      size_ = 0;
      order_ = order;
      builder_.template build<geometryId>(order_);
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
      const Dune::GeometryType geoType = builder_.type();

      std::vector< Field > testBasisVal;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // boundary dofs:
      typedef Dune::QuadratureRule<Field, dimension-1> FaceQuadrature;
      typedef Dune::QuadratureRules<Field, dimension-1> FaceQuadratureRules;

      const auto &refElement = Dune::ReferenceElements< Field, dimension >::general( geoType );

      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        if (!builder_.testFaceBasis(f))
          continue;
        testBasisVal.resize(builder_.testFaceBasis(f)->size());

        const auto &geometry = refElement.template geometry< 1 >( f );
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
    /** \brief evaluate boundary functionals
     *
     *  \param startRow     row of matrix to start
     *  \param mVal         value of the testBasis at a quadrature point on a face
     *  \param rtVal        value of the RaviartThomasBasis at a quadrature point on a face
     *  \param normal       the normal of the face
     *  \param weight       quadrature weight
     *  \param matrix       result gets written into matrix starting with row: row
     */
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
    /** \brief evaluate interior functionals
     *
     *  \param startRow     row of matrix to start
     *  \param mVal         value of the testBasis at a quadrature point in the interior of the ReferenceElement
     *  \param rtVal        value of the RaviartThomasBasis at a quadrature point in the interior of the ReferenceElement
     *  \param weight       quadrature weight
     *  \param matrix       result gets written into matrix starting with row: row
     */
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
    std::size_t order_;
    std::size_t size_;
  };

  template < unsigned int dim, class Field >
  struct RaviartThomasL2InterpolationFactory
  {
    typedef RTL2InterpolationBuilder<dim,Field> Builder;
    typedef const RaviartThomasL2Interpolation<dim,Field> Object;
    typedef std::size_t Key;
    typedef typename std::remove_const<Object>::type NonConstObject;

    template <GeometryType::Id geometryId>
    static Object *create( const Key &key )
    {
      if ( !supports<geometryId>(key) )
        return 0;
      NonConstObject *interpol = new NonConstObject();
      interpol->template build<geometryId>(key);
      return interpol;
    }
    template< GeometryType::Id geometryId >
    static bool supports ( const Key &key )
    {
      return GeometryType(geometryId).isSimplex();
    }
    static void release( Object *object ) { delete object; }
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASSIMPLEX_RAVIARTTHOMASSIMPLEXINTERPOLATION_HH
