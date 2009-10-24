// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASINTERPOLATION_HH
#define DUNE_RAVIARTTHOMASINTERPOLATION_HH
#include <fstream>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/finiteelements/generic/interpolhelper.hh>
#include <dune/finiteelements/common/matrix.hh>
#include <dune/finiteelements/common/localinterpolation.hh>
#include <dune/finiteelements/common/localcoefficients.hh>

#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/quadrature/subquadrature.hh>

#include <dune/finiteelements/generic/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>


namespace Dune
{
  // LocalCoefficientsContainer
  // -------------------
  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactory;
  template < unsigned int dim, class Field >
  struct RaviartThomasL2InterpolationFactory;

  class LocalCoefficientsContainer
    : public LocalCoefficientsInterface< LocalCoefficientsContainer >
  {
    typedef LocalCoefficientsContainer This;
    typedef LocalCoefficientsInterface< This > Base;

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

  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactoryTraits
  {
    static const unsigned int dimension = dim;
    typedef const LocalCoefficientsContainer Object;
    typedef unsigned int Key;
    typedef RaviartThomasCoefficientsFactory<dim> Factory;
  };
  template < unsigned int dim >
  struct RaviartThomasCoefficientsFactory :
    public TopologyFactory< RaviartThomasCoefficientsFactoryTraits<dim> >
  {
    typedef RaviartThomasCoefficientsFactoryTraits<dim> Traits;
    template <class Topology>
    static typename Traits::Object *createObject( const typename Traits::Key &key )
    {
      typedef RaviartThomasL2InterpolationFactory<dim,double> InterpolationFactory;
      if (! supports<Topology>(key) )
        return 0;
      typename InterpolationFactory::Object *interpol
        = InterpolationFactory::template create<Topology>(key);
      typename Traits::Object *localKeys = new typename Traits::Object(*interpol);
      InterpolationFactory::release(interpol);
      return localKeys;
    }
    template< class Topology >
    static bool supports ( const typename Traits::Key &key )
    {
      return GenericGeometry::IsSimplex<Topology>::value;
    }
  };

  // LocalInterpolation
  // -------------------

  // -----------------------------------------
  // RTL2InterpolationBuilder
  // -----------------------------------------
  // L2 Interpolation requires:
  // - for element
  //   - test basis
  // - for each face (dynamic)
  //   - test basis
  //   - normal
  template <unsigned int dim, class Field>
  struct RTL2InterpolationBuilder
  {
    static const unsigned int dimension = dim;
    typedef OrthonormalBasisFactory<dimension,Field> TestBasisFactory;
    typedef typename TestBasisFactory::Object TestBasis;
    typedef OrthonormalBasisFactory<dimension-1,Field> TestFaceBasisFactory;
    typedef typename TestFaceBasisFactory::Object TestFaceBasis;
    typedef FieldVector<Field,dimension> Normal;

    RTL2InterpolationBuilder()
    {}

    unsigned int topologyId() const
    {
      return topologyId_;
    }
    unsigned int order() const
    {
      return order_;
    }
    unsigned int faceSize() const
    {
      return faceSize_;
    }
    TestBasis &testBasis() const
    {
      return *testBasis_;
    }
    TestFaceBasis &testFaceBasis( unsigned int f ) const
    {
      assert( f < faceSize() );
      return *(faceStructure_[f].basis_);
    }
    const Normal &normal( unsigned int f ) const
    {
      return *(faceStructure_[f].normal_);
    }

    template <class Topology>
    void build(unsigned int order)
    {
      order_ = order;
      topologyId_ = Topology::id;
      testBasis_ = TestBasisFactory::template create<Topology>(order-1);
      const unsigned int size = GenericGeometry::Size<Topology,1>::value;
      faceSize_ = size;
      faceStructure_.reserve( faceSize_ );
      ForLoop< Creator<Topology>::template GetCodim,0,size-1>::apply(order, faceStructure_ );
      assert( faceStructure_.size() == faceSize_ );
    }

  private:
    struct FaceStructure
    {
      FaceStructure( TestFaceBasis *tfb,
                     const Normal &n )
        : basis_(tfb), normal_(&n)
      {}
      TestFaceBasis *basis_;
      const Dune::FieldVector<Field,dimension> *normal_;
    };
    template < class Topology >
    struct Creator
    {
      template < int face >
      struct GetCodim
      {
        typedef typename GenericGeometry::SubTopology<Topology,1,face>::type FaceTopology;
        static void apply( const unsigned int order,
                           std::vector<FaceStructure> &faceStructure )
        {
          faceStructure.push_back(
            TestFaceBasisFactory::template create<FaceTopology>(order),
            GenericGeometry::ReferenceElement<Topology,Field>::integrationOuterNormal(face) );
        }
      };
    };

    std::vector<FaceStructure> faceStructure_;
    TestBasis *testBasis_;
    unsigned int topologyId_, order_, faceSize_;
  };

  // A L2 based interpolation for Raviart Thomas
  // --------------------------------------------------
  template< unsigned int dimension, class F>
  class RaviartThomasL2Interpolation
    : public InterpolationHelper<F,dimension>
  {
    typedef RaviartThomasL2Interpolation< dimension, F > This;
    typedef InterpolationHelper<F,dimension> Base;

  public:
    typedef F Field;
    typedef RTL2InterpolationBuilder<dimension,Field> Builder;
    RaviartThomasL2Interpolation( )
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
      order_ = order;
      builder_.template build<Topology>(order_);
      size_ = dimension*builder_.testBasis().size();
      for ( unsigned int f=0; f<builder_.faceSize(); ++f )
        size_ += builder_.testFaceBasis(f).size();
    }

    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      unsigned int row = 0;
      for (unsigned int f=0; f<builder_.faceSize(); ++f)
        for (unsigned int i=0; i<builder_.testFaceBasis(f).size(); ++i,++row)
          keys[row] = LocalKey(f,1,i);
      for (unsigned int i=0; i<builder_.testBasis().size()*dimension; ++i,++row)
        keys[row] = LocalKey(0,0,i);
      assert( row == size() );
    }

  protected:
    template< class Func, class Container, bool type >
    void interpolate ( typename Base::template Helper<Func,Container,type> &func ) const
    {
      const unsigned int topologyId = builder_.topologyId();

      std::vector< Field > testBasisVal;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // boundary dofs:
      typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension-1, Field > FaceQuadratureProvider;
      typedef Dune::GenericGeometry::SubQuadratureProvider< dimension, FaceQuadratureProvider> SubQuadratureProvider;

      for (unsigned int f=0; f<builder_.faceSize(); ++f)
      {
        testBasisVal.resize(builder_.testFaceBasis(f).size());
        const typename SubQuadratureProvider::Quadrature &faceQuad = SubQuadratureProvider::quadrature( topologyId, std::make_pair(f,2*order_+2) );
        const typename SubQuadratureProvider::SubQuadrature &faceSubQuad = SubQuadratureProvider::subQuadrature( topologyId, std::make_pair(f,2*order_+2) );

        const unsigned int quadratureSize = faceQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          builder_.testFaceBasis(f).template evaluate<0>(faceSubQuad.point(qi),testBasisVal);
          fillBnd( row, testBasisVal,
                   func.evaluate(faceQuad.point(qi)),
                   builder_.normal(f), faceQuad.weight(qi),
                   func);
        }
        row += builder_.testFaceBasis(f).size();
      }
      // element dofs
      if (row<size())
      {
        testBasisVal.resize(builder_.testBasis().size());

        typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension, Field > QuadratureProvider;
        const typename QuadratureProvider::Quadrature &elemQuad = QuadratureProvider::template quadrature(topologyId,2*order_+1);
        const unsigned int quadratureSize = elemQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          builder_.testBasis().template evaluate<0>(elemQuad.point(qi),testBasisVal);
          fillInterior( row, testBasisVal,
                        func.evaluate(elemQuad.point(qi)),
                        elemQuad.weight(qi),
                        func );
        }
        row += builder_.testBasis().size()*dimension;
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
    typedef typename remove_const<Object>::type NonConstObject;
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
      return GenericGeometry::IsSimplex<Topology>::value;
    }
  };
}
#endif // DUNE_RAVIARTTHOMASINTERPOLATION_HH
