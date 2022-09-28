// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXINTERPOLATION_HH

#include <fstream>
#include <utility>
#include <numeric>

#include <dune/common/exceptions.hh>

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

  template < unsigned int dim, class Field >
  struct NedelecL2InterpolationFactory;



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



  // NedelecCoefficientsFactory
  // --------------------------------

  template < unsigned int dim >
  struct NedelecCoefficientsFactory
  {
    typedef std::size_t Key;
    typedef const LocalCoefficientsContainer Object;

    template< GeometryType::Id geometryId >
    static Object *create( const Key &key )
    {
      typedef NedelecL2InterpolationFactory< dim, double > InterpolationFactory;
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
      GeometryType gt = geometryId;
      return gt.isTriangle() || gt.isTetrahedron() ;
    }
    static void release( Object *object ) { delete object; }
  };



  // NedelecL2InterpolationBuilder
  // ------------------------

  // L2 Interpolation requires:
  // - for element
  //   - test basis
  // - for each face (dynamic)
  //   - test basis
  //   - tangents
  // - for each edge (dynamic)
  //   - test basis
  //   - tangent
  template< unsigned int dim, class Field >
  struct NedelecL2InterpolationBuilder
  {
    static const unsigned int dimension = dim;

    // for the dofs associated to the element
    typedef OrthonormalBasisFactory< dimension, Field > TestBasisFactory;
    typedef typename TestBasisFactory::Object TestBasis;

    // for the dofs associated to the faces
    typedef OrthonormalBasisFactory< dimension-1, Field > TestFaceBasisFactory;
    typedef typename TestFaceBasisFactory::Object TestFaceBasis;

    // for the dofs associated to the edges
    typedef OrthonormalBasisFactory< 1, Field > TestEdgeBasisFactory;
    typedef typename TestEdgeBasisFactory::Object TestEdgeBasis;

    // the tangent of the edges
    typedef FieldVector< Field, dimension > Tangent;

    // the normal and the tangents of the faces
    typedef FieldVector< Field, dimension > Normal;
    typedef std::array<FieldVector< Field, dimension >,dim-1> FaceTangents;

    NedelecL2InterpolationBuilder () = default;

    NedelecL2InterpolationBuilder ( const NedelecL2InterpolationBuilder & ) = delete;
    NedelecL2InterpolationBuilder ( NedelecL2InterpolationBuilder && ) = delete;

    ~NedelecL2InterpolationBuilder ()
    {
      TestBasisFactory::release( testBasis_ );
      for( FaceStructure &f : faceStructure_ )
        TestFaceBasisFactory::release( f.basis_ );
      for( EdgeStructure& e : edgeStructure_ )
        TestEdgeBasisFactory::release( e.basis_ );
    }

    unsigned int topologyId () const
    {
      return geometry_.id();
    }

    GeometryType type () const
    {
      return geometry_;
    }

    std::size_t order () const
    {
      return order_;
    }

    // number of faces
    unsigned int faceSize () const
    {
      return numberOfFaces_;
    }

    // number of edges
    unsigned int edgeSize () const
    {
      return numberOfEdges_;
    }

    // basis associated to the element
    TestBasis *testBasis () const
    {
      return testBasis_;
    }

    // basis associated to face f
    TestFaceBasis *testFaceBasis ( unsigned int f ) const
    {
      assert( f < faceSize() );
      return faceStructure_[ f ].basis_;
    }

    // basis associated to edge e
    TestEdgeBasis *testEdgeBasis ( unsigned int e ) const
    {
      assert( e < edgeSize() );
      return edgeStructure_[ e ].basis_;
    }

    const Tangent& edgeTangent ( unsigned int e ) const
    {
      assert( e < edgeSize() );
      return edgeStructure_[ e ].tangent_;
    }

    const FaceTangents& faceTangents ( unsigned int f ) const
    {
      assert( f < faceSize() );
      return faceStructure_[ f ].faceTangents_;
    }

    const Normal &normal ( unsigned int f ) const
    {
      assert( f < faceSize() );
      return faceStructure_[ f ].normal_;
    }

    template< GeometryType::Id geometryId >
    void build ( std::size_t order )
    {
      constexpr GeometryType geometry = geometryId;
      order_ = order;
      geometry_ = geometry;

      /*
       * The Nedelec parameter begins at 1.
       * This is the numbering used by J.C. Nedelec himself.
       * See "Mixed Finite Elements in \R^3" published in 1980.
       *
       * This construction is based on the construction of Raviart-Thomas elements.
       * There the numbering starts at 0.
       * Because of this we reduce the order internally by 1.
       */
      order--;

      // if dimension == 2: order-1 on element
      // if dimension == 3: order-2 on element
      int requiredOrder =  static_cast<int>(dimension==3);
      testBasis_ = (order > requiredOrder ? TestBasisFactory::template create< geometry >( order-1-requiredOrder ) : nullptr);

      const auto &refElement = ReferenceElements< Field, dimension >::general( type() );

      numberOfFaces_ = refElement.size( 1 );
      faceStructure_.reserve( numberOfFaces_ );

      // compute the basis, tangents and normals of each face
      for (std::size_t i=0; i<numberOfFaces_; i++)
      {
        FieldVector<Field,dimension> zero(0);
        FaceTangents faceTangents;
        faceTangents.fill(zero);

        // use the first dim-1 vertices of a face to compute the tangents
        auto vertices = refElement.subEntities(i,1,dim).begin();
        auto vertex1 = *vertices;
        for(int j=1; j<dim;j++)
        {
          auto vertex2 = vertices[j];

          faceTangents[j-1] = refElement.position(vertex2,dim) - refElement.position(vertex1,dim);

          // By default, edges point from the vertex with the smaller index
          // to the vertex with the larger index.
          if (vertex1>vertex2)
            faceTangents[j-1] *=-1;

          vertex1 = vertex2;
        }

        /* For simplices or cubes of arbitrary dimension you could just use
         *
         * ```
         * GeometryType faceGeometry = Impl::getBase(geometry_);
         * TestFaceBasis *faceBasis = ( dim == 3 && order > 0 ? TestFaceBasisFactory::template create< faceGeometry >( order-1 ) : nullptr);
         * ```
         *
         * For i.e. Prisms and Pyramids in 3d this does not work because they contain squares and triangles as faces.
         * And depending on the dynamic face index a different face geometry is needed.
         *
         */
        TestFaceBasis *faceBasis = ( dim == 3 && order > 0 ? Impl::IfGeometryType< CreateFaceBasis, dimension-1 >::apply( refElement.type( i, 1 ), order-1 ) : nullptr);
        faceStructure_.emplace_back( faceBasis, refElement.integrationOuterNormal(i), faceTangents );
      }
      assert( faceStructure_.size() == numberOfFaces_ );

      numberOfEdges_ = refElement.size( dim-1 );
      edgeStructure_.reserve( numberOfEdges_ );

      // compute the basis and tangent of each edge
      for (std::size_t i=0; i<numberOfEdges_; i++)
      {
        auto vertexIterator = refElement.subEntities(i,dim-1,dim).begin();
        auto v0 = *vertexIterator;
        auto v1 = *(++vertexIterator);

        // By default, edges point from the vertex with the smaller index
        // to the vertex with the larger index.
        if (v0>v1)
          std::swap(v0,v1);
        auto tangent = std::move(refElement.position(v1,dim) - refElement.position(v0,dim));

        TestEdgeBasis *edgeBasis = Impl::IfGeometryType< CreateEdgeBasis, 1 >::apply( refElement.type( i, dim-1 ), order );
        edgeStructure_.emplace_back( edgeBasis, tangent );
      }
      assert( edgeStructure_.size() == numberOfEdges_ );
    }

  private:

    // helper struct for edges
    struct EdgeStructure
    {
      EdgeStructure( TestEdgeBasis *teb, const Tangent &t )
        : basis_( teb ), tangent_( t )
      {}

      TestEdgeBasis *basis_;
      const Dune::FieldVector< Field, dimension > tangent_;
    };

    template< GeometryType::Id edgeGeometryId >
    struct CreateEdgeBasis
    {
      static TestEdgeBasis *apply ( std::size_t order ) { return TestEdgeBasisFactory::template create< edgeGeometryId >( order ); }
    };

    // helper struct for faces
    struct FaceStructure
    {
      FaceStructure( TestFaceBasis *tfb, const Normal& normal, const FaceTangents& faceTangents )
        : basis_( tfb ), normal_(normal), faceTangents_( faceTangents )
      {}

      TestFaceBasis *basis_;
      const Dune::FieldVector< Field, dimension > normal_;
      const FaceTangents faceTangents_;
    };

    template< GeometryType::Id faceGeometryId >
    struct CreateFaceBasis
    {
      static TestFaceBasis *apply ( std::size_t order ) { return TestFaceBasisFactory::template create< faceGeometryId >( order ); }
    };

    TestBasis *testBasis_ = nullptr;
    std::vector< FaceStructure > faceStructure_;
    unsigned int numberOfFaces_;
    std::vector< EdgeStructure > edgeStructure_;
    unsigned int numberOfEdges_;
    GeometryType geometry_;
    std::size_t order_;
  };



  // NedelecL2Interpolation
  // ----------------------------

  /**
   * \class NedelecL2Interpolation
   * \brief An L2-based interpolation for Nedelec
   *
   **/
  template< unsigned int dimension, class F>
  class NedelecL2Interpolation
    : public InterpolationHelper< F ,dimension >
  {
    typedef NedelecL2Interpolation< dimension, F > This;
    typedef InterpolationHelper<F,dimension> Base;

  public:
    typedef F Field;
    typedef NedelecL2InterpolationBuilder<dimension,Field> Builder;
    typedef typename Builder::FaceTangents FaceTangents;

    NedelecL2Interpolation()
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
          size_ += (dimension-1)*builder_.testFaceBasis(f)->size();

      for ( unsigned int e=0; e<builder_.edgeSize(); ++e )
        if (builder_.testEdgeBasis(e))
          size_ += builder_.testEdgeBasis(e)->size();
    }

    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      unsigned int row = 0;
      for (unsigned int e=0; e<builder_.edgeSize(); ++e)
      {
        if (builder_.edgeSize())
          for (unsigned int i=0; i<builder_.testEdgeBasis(e)->size(); ++i,++row)
            keys[row] = LocalKey(e,dimension-1,i);
      }
      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        if (builder_.faceSize())
          for (unsigned int i=0; i<builder_.testFaceBasis(f)->size()*(dimension-1); ++i,++row)
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

      std::vector<Field> testBasisVal;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // edge dofs:
      typedef Dune::QuadratureRule<Field, 1> EdgeQuadrature;
      typedef Dune::QuadratureRules<Field, 1> EdgeQuadratureRules;

      const auto &refElement = Dune::ReferenceElements< Field, dimension >::general( geoType );

      for (unsigned int e=0; e<builder_.edgeSize(); ++e)
      {
        if (!builder_.testEdgeBasis(e))
          continue;
        testBasisVal.resize(builder_.testEdgeBasis(e)->size());

        const auto &geometry = refElement.template geometry< dimension-1 >( e );
        const Dune::GeometryType subGeoType( geometry.type().id(), 1 );
        const EdgeQuadrature &edgeQuad = EdgeQuadratureRules::rule( subGeoType, 2*order_+2 );

        const unsigned int quadratureSize = edgeQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          if (dimension>1)
            builder_.testEdgeBasis(e)->template evaluate<0>(edgeQuad[qi].position(),testBasisVal);
          else
            testBasisVal[0] = 1.;
          computeEdgeDofs(row,
                          testBasisVal,
                          func.evaluate( geometry.global( edgeQuad[qi].position() ) ),
                          builder_.edgeTangent(e),
                          edgeQuad[qi].weight(),
                          func);
        }

        row += builder_.testEdgeBasis(e)->size();
      }

      // face dofs:
      typedef Dune::QuadratureRule<Field, dimension-1> FaceQuadrature;
      typedef Dune::QuadratureRules<Field, dimension-1> FaceQuadratureRules;

      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        if (builder_.testFaceBasis(f))
        {
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

            computeFaceDofs( row,
                             testBasisVal,
                             func.evaluate( geometry.global( faceQuad[qi].position() ) ),
                             builder_.faceTangents(f),
                             builder_.normal(f),
                             faceQuad[qi].weight(),
                             func);
          }

          row += builder_.testFaceBasis(f)->size()*(dimension-1);
        }
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
          computeInteriorDofs(row,
                              testBasisVal,
                              func.evaluate(elemQuad[qi].position()),
                              elemQuad[qi].weight(),
                              func );
        }

        row += builder_.testBasis()->size()*dimension;
      }
      assert(row==size());
    }

  private:
    /** /brief evaluate functionals associated to the edges
     *
     *  \param startRow    row of matrix to start
     *  \param mVal         value of the testBasis at a quadrature point on an edge
     *  \param nedVal       value of the nedelecBasis at a quadrature point on an edge
     *  \param tangent      the tangent of the edge
     *  \param weight       quadrature weight
     *  \param matrix       result gets written into matrix starting with row: row
     */
    template <class MVal, class NedVal,class Matrix>
    void computeEdgeDofs (unsigned int startRow,
                          const MVal &mVal,
                          const NedVal &nedVal,
                          const FieldVector<Field,dimension> &tangent,
                          const Field &weight,
                          Matrix &matrix) const
    {
      const unsigned int endRow = startRow+mVal.size();
      typename NedVal::const_iterator nedIter = nedVal.begin();
      for ( unsigned int col = 0; col < nedVal.size() ; ++nedIter,++col)
      {
        Field cFactor = (*nedIter)*tangent;
        typename MVal::const_iterator mIter = mVal.begin();
        for (unsigned int row = startRow; row!=endRow; ++mIter, ++row )
          matrix.add(row,col, (weight*cFactor)*(*mIter) );

        assert( mIter == mVal.end() );
      }
    }

    /** /brief evaluate functionals associated to the faces
     *
     *  \param startRow      row of matrix to start
     *  \param mVal           value of the testBasis at a quadrature point on a face
     *  \param nedVal         value of the nedelecBasis at a quadrature point on a face
     *  \param faceTangents   the tangents of the face
     *  \param normal         the normal of the face
     *  \param weight         quadrature weight
     *  \param matrix         result gets written into matrix starting with row: row
     */
    template <class MVal, class NedVal,class Matrix>
    void computeFaceDofs (unsigned int startRow,
                          const MVal &mVal,
                          const NedVal &nedVal,
                          const FaceTangents& faceTangents,
                          const FieldVector<Field,dimension> &normal,
                          const Field &weight,
                          Matrix &matrix) const
    {
      const unsigned int endRow = startRow+mVal.size()*(dimension-1);
      typename NedVal::const_iterator nedIter = nedVal.begin();
      for ( unsigned int col = 0; col < nedVal.size() ; ++nedIter,++col)
      {
        auto const& u=*nedIter;
        auto const& n=normal;
        FieldVector<Field,dimension> nedTimesNormal = { u[1]*n[2]-u[2]*n[1],
                                                        u[2]*n[0]-u[0]*n[2],
                                                        u[0]*n[1]-u[1]*n[0]};
        typename MVal::const_iterator mIter = mVal.begin();
        for (unsigned int row = startRow; row!=endRow; ++mIter)
        {
          for(int i=0; i<dimension-1;i++)
          {
            auto test = *mIter*faceTangents[i];
            matrix.add(row+i,col, weight*(nedTimesNormal*test) );
          }
          row += dimension-1;
        }

        assert( mIter == mVal.end() );
      }
    }

    /** /brief evaluate functionals associated to the interior
     *
     *  \param startRow       row of matrix to start
     *  \param mVal           value of the testBasis at a quadrature point in the interior of the ReferenceElement
     *  \param nedVal         value of the nedelecBasis at a quadrature point in the interior of the ReferenceElement
     *  \param weight         quadrature weight
     *  \param matrix         result gets written into matrix starting with row: row
     */
    template <class MVal, class NedVal,class Matrix>
    void computeInteriorDofs (unsigned int startRow,
                              const MVal &mVal,
                              const NedVal &nedVal,
                              Field weight,
                              Matrix &matrix) const
    {
      const unsigned int endRow = startRow+mVal.size()*dimension;
      typename NedVal::const_iterator nedIter = nedVal.begin();
      for ( unsigned int col = 0; col < nedVal.size() ; ++nedIter,++col)
      {
        typename MVal::const_iterator mIter = mVal.begin();
        for (unsigned int row = startRow; row!=endRow; ++mIter,row+=dimension )
          for (unsigned int i=0; i<dimension; ++i)
            matrix.add(row+i,col, (weight*(*mIter))*(*nedIter)[i] );

        assert( mIter == mVal.end() );
      }
    }

  public:
    Builder builder_;
    std::size_t order_;
    std::size_t size_;
  };

  template < unsigned int dim, class Field >
  struct NedelecL2InterpolationFactory
  {
    typedef NedelecL2InterpolationBuilder<dim,Field> Builder;
    typedef const NedelecL2Interpolation<dim,Field> Object;
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

    template <GeometryType::Id geometryId>
    static bool supports( const Key &key )
    {
      GeometryType gt = geometryId;
      return gt.isTriangle() || gt.isTetrahedron() ;
    }
    static void release( Object *object ) { delete object; }
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXINTERPOLATION_HH
