// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASBASIS_HH
#define DUNE_RAVIARTTHOMASBASIS_HH
#include <fstream>
#include <utility>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>

#include <dune/finiteelements/common/localcoefficients.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/lobattopoints.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/quadrature/subquadrature.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>
#include <dune/finiteelements/lagrangebasis/lobattopoints.hh>

namespace Dune
{

  // LocalCoefficientsContainer
  // -------------------

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

  // A small helper class to avoid having to
  // write the interpolation twice (once for function
  // and once for a basis)
  template< class F, unsigned int dimension >
  struct RaviartThomasInterpolation
  {
    template <class Func,class Container, bool type>
    struct Helper;
  };
  template <class F,unsigned int d>
  template <class Func,class Vector>
  struct RaviartThomasInterpolation<F,d>::Helper<Func,Vector,true>
  // Func is of Function type
  {
    typedef std::vector< Dune::FieldVector<F,d> > Result;
    Helper(const Func & func, Vector &vec)
      : func_(func),
        vec_(vec),
        tmp_(1)
    {}
    const typename Vector::value_type &operator()(unsigned int row,unsigned int col)
    {
      return vec_[row];
    }
    template <class Fy>
    void set(unsigned int row,unsigned int col,
             const Fy &val)
    {
      assert(-100<val && val<100);
      assert(col==0);
      assert(row<vec_.size());
      field_cast( val, vec_[row] );
    }
    template <class Fy>
    void add(unsigned int row,unsigned int col,
             const Fy &val)
    {
      assert(-100<val && val<100);
      assert(col==0);
      assert(row<vec_.size());
      vec_[row] += field_cast<typename Vector::value_type>(val);
    }
    template <class DomainVector>
    const Result &evaluate(const DomainVector &x) const
    {
      field_cast(func_( x ), tmp_[0] );
      return tmp_;
    }
    unsigned int size() const
    {
      return 1;
    }
    const Func &func_;
    Vector &vec_;
    mutable Result tmp_;
  };
  template <class F,unsigned int d>
  template <class Basis,class Matrix>
  struct RaviartThomasInterpolation<F,d>::Helper<Basis,Matrix,false>
  // Func is of Basis type
  {
    typedef std::vector< Dune::FieldVector<F,d> > Result;
    Helper(const Basis & basis, Matrix &matrix)
      : basis_(basis),
        matrix_(matrix),
        tmp_(basis.size()) {}
    const F &operator()(unsigned int row,unsigned int col) const
    {
      return matrix_(row,col);
    }
    F &operator()(unsigned int row,unsigned int col)
    {
      return matrix_(row,col);
    }
    template <class Fy>
    void set(unsigned int row,unsigned int col,
             const Fy &val)
    {
      assert(col<matrix_.cols());
      assert(row<matrix_.rows());
      field_cast(val,matrix_(row,col));
    }
    template <class Fy>
    void add(unsigned int row,unsigned int col,
             const Fy &val)
    {
      assert(col<matrix_.cols());
      assert(row<matrix_.rows());
      matrix_(row,col) += val;
    }
    template <class DomainVector>
    const Result &evaluate(const DomainVector &x) const
    {
      basis_.template evaluate<0>(x,tmp_);
      return tmp_;
    }
    unsigned int size() const
    {
      return basis_.size();
    }

    const Basis &basis_;
    Matrix &matrix_;
    mutable Result tmp_;
  };

  // -----------------------------------------
  // Basis
  // -----------------------------------------
  template <class Field,unsigned int dim>
  struct RaviartThomasInitialBasis
  {
    static const unsigned int dimension = dim;
    // typedef MonomialBasisProvider<dimension,Field> TestBasisProvider;
    // typedef MonomialBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    typedef OrthonormalBasisProvider<dimension,Field> TestBasisProvider;
    typedef OrthonormalBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef LagrangeBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LagrangeBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef LobattoBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LobattoBasisProvider<dimension-1,Field> TestFaceBasisProvider;

    typedef typename TestBasisProvider::Basis TestBasis;
    typedef typename TestFaceBasisProvider::Basis TestFaceBasis;

    struct FaceStructure
    {
      template < class Topology >
      struct Get
      {
        static const TestBasis &testBasis(unsigned int order)
        {
          return TestBasisProvider::template basis<Topology>(order-1);
        }
        template < int face >
        struct GetCodim
        {
          typedef typename GenericGeometry::SubTopology<Topology,1,face>::type FaceTopology;
          static void apply( const unsigned int order,
                             std::vector<FaceStructure*> &faceStructure )
          {
            faceStructure.push_back( new FaceStructure (
                                       TestFaceBasisProvider::template basis<FaceTopology>(order),
                                       GenericGeometry::ReferenceElement<Topology,Field>::integrationOuterNormal(face) ) );
          }
        };
      };
      FaceStructure( const TestFaceBasis &tfb, const Dune::FieldVector<Field,dimension> &n )
        : basis(tfb), normal(n)
      {}
      const TestFaceBasis &basis;
      const Dune::FieldVector<Field,dimension> &normal;
    };

  };

  // A L2 based interpolation for Raviart Thomas
  // --------------------------------------------------
  template< class F, unsigned int dimension>
  class RaviartThomasL2Interpolation
    : public RaviartThomasInterpolation<F,dimension>
  {
    typedef RaviartThomasL2Interpolation< F, dimension > This;
    typedef RaviartThomasInterpolation<F,dimension> Base;

  public:
    typedef F Field;
    typedef MonomialBasisProvider<dimension,Field> TestMBasisProvider;
    typedef MonomialBasisProvider<dimension-1,Field> TestFaceMBasisProvider;
    typedef SparseCoeffMatrix< Field, 1 > CoefficientMatrix;
    typedef StandardEvaluator< typename TestMBasisProvider::Basis > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > TestBasis;
    typedef StandardEvaluator< typename TestFaceMBasisProvider::Basis > FaceEvaluator;
    typedef PolynomialBasis< FaceEvaluator, CoefficientMatrix > TestFaceBasis;
    typedef typename RaviartThomasInitialBasis<Field,dimension> :: FaceStructure FaceStructure;
    RaviartThomasL2Interpolation
      ( const unsigned int topologyId,
      const unsigned int order,
      const TestBasis &mBasis,
      const std::vector< FaceStructure* > &faceStructure )
      : order_(order),
        faceStructure_(faceStructure),
        topologyId_(topologyId),
        mBasis_(mBasis),
        size_( 0 )
    {
      size_ = dimension*mBasis_.size();
      for ( unsigned int f=0; f<faceStructure_.size(); ++f )
        size_ += faceStructure_[f]->basis.size();
    }

    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return size_;
    }

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

    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      unsigned int row = 0;
      for (unsigned int f=0; f<dimension+1; ++f)
        for (unsigned int i=0; i<faceStructure_[f]->basis.size(); ++i,++row)
          keys[row] = LocalKey(f,1,i);
      for (unsigned int i=0; i<mBasis_.size()*dimension; ++i,++row)
        keys[row] = LocalKey(0,0,i);
      assert( row == size() );
    }

  protected:
    template< class Func, class Container, bool type >
    void interpolate ( typename Base::template Helper<Func,Container,type> &func ) const
    {
      std::vector< Field > testBasisVal;
      std::vector< Field > testBasisInt;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // boundary dofs:
      typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension-1, Field > FaceQuadratureProvider;
      typedef Dune::GenericGeometry::SubQuadratureProvider< dimension, FaceQuadratureProvider> SubQuadratureProvider;

      for (unsigned int f=0; f<faceStructure_.size(); ++f)
      {
        testBasisVal.resize(faceStructure_[f]->basis.size());
        // testBasisInt.resize(mFaceBasis_.size());
        // mFaceBasis_.integral(testBasisInt);
        const typename SubQuadratureProvider::Quadrature &faceQuad = SubQuadratureProvider::quadrature( topologyId_, std::make_pair(f,2*order_+2) );
        const typename SubQuadratureProvider::SubQuadrature &faceSubQuad = SubQuadratureProvider::subQuadrature( topologyId_, std::make_pair(f,2*order_+2) );

        const unsigned int quadratureSize = faceQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          faceStructure_[f]->basis.template evaluate<0>(faceSubQuad.point(qi),testBasisVal);
          // for (int i=0;i<testBasisVal.size();++i)
          //  testBasisVal[i] /= testBasisInt[i];
          fillBnd( row, testBasisVal,
                   func.evaluate(faceQuad.point(qi)),
                   faceStructure_[f]->normal, faceQuad.weight(qi),
                   func);
        }
        row += faceStructure_[f]->basis.size();
      }
      // element dofs
      if (row<size())
      {
        testBasisVal.resize(mBasis_.size());
        testBasisInt.resize(mBasis_.size());

        typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension, Field > QuadratureProvider;
        const typename QuadratureProvider::Quadrature &elemQuad = QuadratureProvider::template quadrature(topologyId_,2*order_+1);
        const unsigned int quadratureSize = elemQuad.size();
        // mBasis_.integral(testBasisInt);
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          mBasis_.template evaluate<0>(elemQuad.point(qi),testBasisVal);
          // for (int i=0;i<testBasisVal.size();++i)
          //   testBasisVal[i] /= testBasisInt[i];
          fillInterior( row, testBasisVal,
                        func.evaluate(elemQuad.point(qi)),
                        elemQuad.weight(qi),
                        func );
        }
        row += mBasis_.size()*dimension;
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

    unsigned int order_;
    const std::vector< FaceStructure* > faceStructure_;
    const unsigned int topologyId_;
    const TestBasis &mBasis_;
    unsigned int size_;
  };

  template <class Topology, class Field>
  struct RTVecMatrix
  {
    enum {dim = Topology::dimension};
    typedef MultiIndex<dim> MI;
    typedef MonomialBasis<Topology,MI> MIBasis;
    RTVecMatrix(unsigned int order)
    {
      MIBasis basis(order+1);
      FieldVector< MI, dim > x;
      for( int i = 0; i < dim; ++i )
        x[ i ].set( i, 1 );
      std::vector< MI > val( basis.size() );
      basis.evaluate( x, val );

      col_ = basis.size();
      unsigned int notHomogen = 0;
      if (order>0)
        notHomogen = basis.sizes()[order-1];
      unsigned int homogen = basis.sizes()[order]-notHomogen;
      row_ = (notHomogen*dim+homogen*(dim+1))*dim;
      row1_ = basis.sizes()[order]*dim*dim;
      mat_ = new Field*[row_];
      int row = 0;
      for (unsigned int i=0; i<notHomogen+homogen; ++i)
      {
        for (unsigned int r=0; r<dim; ++r)
        {
          for (unsigned int rr=0; rr<dim; ++rr)
          {
            mat_[row] = new Field[col_];
            for (unsigned int j=0; j<col_; ++j)
            {
              mat_[row][j] = 0.;
            }
            if (r==rr)
              mat_[row][i] = 1.;
            ++row;
          }
        }
      }
      for (unsigned int i=0; i<homogen; ++i)
      {
        for (unsigned int r=0; r<dim; ++r)
        {
          mat_[row] = new Field[col_];
          for (unsigned int j=0; j<col_; ++j)
          {
            mat_[row][j] = 0.;
          }
          unsigned int w;
          MI xval = val[notHomogen+i];
          xval *= x[r];
          for (w=homogen+notHomogen; w<val.size(); ++w)
          {
            if (val[w] == xval)
            {
              mat_[row][w] = 1.;
              break;
            }
          }
          assert(w<val.size());
          ++row;
        }
      }
    }
    ~RTVecMatrix()
    {
      for (unsigned int i=0; i<rowSize(); ++i) {
        delete [] mat_[i];
      }
      delete [] mat_;
    }
    unsigned int colSize(int row) const {
      return col_;
    }
    unsigned int rowSize() const {
      return row_;
    }
    const Field operator() ( int r, int c ) const
    {
      assert(r<(int)row_ && c<(int)col_);
      return mat_[r][c];
    }
    unsigned int row_,col_,row1_;
    Field **mat_;
  };

  template <class Topology, class RTInterpolation>
  struct RaviartThomasMatrix {
    typedef typename RTInterpolation::Field Field;
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< Field > mat_t;
    typedef MonomialBasis<Topology,Field> MBasis;
    RaviartThomasMatrix(const RTInterpolation &interpolation) :
      order_(interpolation.order()),
      interpolation_(interpolation),
      vecMatrix_(order_)
    {
      MBasis basis(order_+1);
      typedef StandardEvaluator<MBasis> EvalMBasis;
      typedef PolynomialBasisWithMatrix<EvalMBasis,SparseCoeffMatrix<Field,dim> > TMBasis;
      TMBasis tmBasis(basis);
      tmBasis.fill(vecMatrix_);
      interpolation_.interpolate( tmBasis , matrix_ );
      matrix_.invert();
    }
    unsigned int colSize(int row) const {
      return vecMatrix_.colSize(row);
    }
    unsigned int rowSize() const {
      return vecMatrix_.rowSize();
    }
    const Field operator() ( int r, int c ) const
    {
      int rmod = r%dim;
      int rdiv = r/dim;
      Field ret = 0;
      for (unsigned int k=0; k<matrix_.cols(); ++k) {
        ret += matrix_(k,rdiv)*vecMatrix_(k*dim+rmod,c);
      }
      return ret;
    }
    int order_;
    const RTInterpolation &interpolation_;
    RTVecMatrix<Topology,Field> vecMatrix_;
    mat_t matrix_;
  };


  template< int dim, class SF, class CF >
  struct RaviartThomasBasisCreator
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasis>,SparseCoeffMatrix<StorageField,dim> > Basis;
    typedef LocalCoefficientsContainer LocalCoefficients;

    typedef unsigned int Key;
    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    typedef RaviartThomasL2Interpolation< ComputeField, dimension > LocalInterpolation;

    template< class Topology >
    static const LocalInterpolation &localInterpolation ( const Key &order )
    {
      typedef RaviartThomasInitialBasis<ComputeField,dimension> InitialBasis;
      typedef typename InitialBasis::FaceStructure::template Get<Topology> BasisGet;
      const unsigned int size = GenericGeometry::Size<Topology,1>::value;
      std::vector< typename InitialBasis::FaceStructure* > faceStructure;
      GenericGeometry::ForLoop< BasisGet::template GetCodim,0,size-1>::apply(order, faceStructure );

      const typename InitialBasis::TestBasis &testBasis( BasisGet::testBasis(order) );
      LocalInterpolation *interpolation = new LocalInterpolation( Topology::id,order,testBasis,faceStructure );
      return *interpolation;
    }
    template <class Topology>
    static Basis &basis(unsigned int order)
    {
      const MBasis &_mBasis = MonomialBasisProvider<dimension,StorageField>::template basis<SimplexTopology>(order+1);
      Basis *basis = new Basis(_mBasis);
      const LocalInterpolation &interpolation = localInterpolation<Topology>(order);
      RaviartThomasMatrix<Topology,LocalInterpolation> matrix(interpolation);
      release(interpolation);
      basis->fill(matrix);
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,dimension> > BasisMI;
        const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<SimplexTopology>(order+1);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "rt_" << Topology::name() << "_p" << order;
        std::ofstream out(name.str().c_str());
        basisPrint<0>(out,basisMI);
      }
      return *basis;
    }
    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      const LocalInterpolation &interpolation = localInterpolation<Topology>(order);
      LocalCoefficientsContainer *localKeys = new LocalCoefficientsContainer(interpolation);
      release(interpolation);
      return *localKeys;
    }

    static void release ( const LocalInterpolation &localInterpolation )
    {
      delete &localInterpolation;
    }
    static void release ( const Basis &basis )
    {
      delete &basis;
    }
    static void release ( const LocalCoefficients &localCoefficients )
    {
      delete &localCoefficients;
    }
    template< class Topology >
    static bool supports ( const Key &key )
    {
      return GenericGeometry::IsSimplex<Topology>::value;
    }
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct RaviartThomasBasisProvider
    : public BasisProvider<RaviartThomasBasisCreator<dim,SF,CF> >
  {};

}
#endif // DUNE_RAVIARTTHOMASBASIS_HH
