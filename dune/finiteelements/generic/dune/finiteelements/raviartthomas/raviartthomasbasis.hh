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
  // AllBasisCreator
  // ---------------
  template <unsigned int dimension,class Field>
  struct ScalarBasisProvider
  {
    enum Type {orthonormal,lagrange,lobatto};
    typedef OrthonormalBasisProvider<dimension,Field> ONBProvider;
    typedef LagrangeBasisProvider<dimension,Field> LagrangeProvider;
    typedef LobattoBasisProvider<dimension,Field> LobattoProvider;

    typedef MonomialBasisProvider<dimension,Field> MBasisProvider;
    typedef SparseCoeffMatrix< Field, 1 > CoefficientMatrix;
    typedef StandardEvaluator< typename MBasisProvider::Basis > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > Basis;
    template <class Topology>
    static const Basis& basis(Type type,unsigned int order)
    {
      switch (type)
      {
      case orthonormal :
        return ONBProvider::template basis<Topology>(order);
      case lagrange :
        return LagrangeProvider::template basis<Topology>(order);
      case lobatto :
        return LobattoProvider::template basis<Topology>(order);
      }
    }
    static void release(Type type,const Basis& basis)
    {
      switch (type)
      {
      case orthonormal :
        ONBProvider::release( static_cast<const typename ONBProvider::Basis&>(basis) );
        break;
      case lagrange :
        LagrangeProvider::release( static_cast<const typename LagrangeProvider::Basis&>(basis) );
        break;
      case lobatto :
        LobattoProvider::release( static_cast<const typename LobattoProvider::Basis&>(basis) );
        break;
      }
    }
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

  // A lagrange based interpolation for Raviart Thomas
  // --------------------------------------------------
  template< class F, class PointsSetCreator >
  class RaviartThomasLagrangeInterpolation
    : public RaviartThomasInterpolation<F,PointsSetCreator::LagrangePoints::dimension>
  {
    typedef typename PointsSetCreator::LagrangePoints LagrangePoints;
    static const unsigned int dimension = LagrangePoints::dimension;

    typedef RaviartThomasLagrangeInterpolation< F, PointsSetCreator > This;
    typedef RaviartThomasInterpolation<F,dimension> Base;

  public:
    typedef F Field;

  public:
    RaviartThomasLagrangeInterpolation
      ( const unsigned int order,
      const LagrangePoints &points,
      const FieldMatrix<Field,dimension+1,dimension> &normal )
      : order_(order),
        lagrangePoints_ ( points ),
        normal_(normal),
        size_(0)
    {
      typedef typename LagrangePoints::iterator Iterator;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
      {
        if (it->localKey().codim()==1)
          size_ += 1;
        else if (it->localKey().codim()==0)
          size_ += dimension;
      }
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
      typename Base::template Helper<Function,std::vector<Fy>,true> func( function, coefficients );
      interpolate(func);
    }
    template< class Basis, class Matrix >
    void interpolate ( const Basis &basis, Matrix &matrix ) const
    {
      matrix.resize( size(), basis.size() );
      typename Base::template Helper<Basis,Matrix,false> func( basis, matrix );
      interpolate(func);
    }

    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      typedef typename LagrangePoints::iterator Iterator;
      unsigned int row = 0;
      unsigned int interior = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
      {
        if (it->localKey().codim()==1)
        {
          keys[row] = it->localKey();
          ++row;
        }
        else if (it->localKey().codim()==0)
        {
          for (int i=0; i<dimension; ++i,++row,++interior)
            keys[row] = LocalKey(0,0,interior);
        }
      }
      assert( row == size() );
    }

  protected:
    template< class Func, class Container, bool type >
    void interpolate ( typename Base::template Helper<Func,Container,type> &func ) const
    {
      typedef typename LagrangePoints::iterator Iterator;
      unsigned int row = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
      {
        if (it->localKey().codim()>1)
          continue;
        if (it->localKey().codim()==1)
          fillBnd ( row, normal_[it->localKey().subEntity()],
                    func.evaluate( it->point() ),
                    func );
        else
          fillInterior ( row,
                         func.evaluate( it->point() ),
                         func );
      }
    }

  private:
    /** /brief evaluate boundary functionals **/
    template <class Val,class Result>
    void fillBnd (unsigned int &row,
                  const FieldVector<Field,dimension> &normal,
                  const Val &val,
                  Result &res) const
    {
      typename Val::const_iterator iter = val.begin();
      for ( unsigned int col = 0; col < val.size() ; ++iter,++col)
        res.set(row,col, (*iter)*normal );
      ++row;
    }
    template <class Val,class Result>
    void fillInterior (unsigned int &row,
                       const Val &val,
                       Result &res) const
    {
      typename Val::const_iterator iter = val.begin();
      for ( unsigned int col = 0; col < val.size() ; ++iter,++col)
        for (unsigned int d=0; d<dimension; ++d)
          res.set(row+d,col, (*iter)[d] );
      row+=dimension;
    }

    unsigned int order_;
    FieldMatrix<Field,dimension+1,dimension> normal_;
    const LagrangePoints &lagrangePoints_;
    unsigned int size_;
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
    typedef typename GenericGeometry::SimplexTopology<dimension>::type Topology;
    typedef typename GenericGeometry::SimplexTopology<dimension-1>::type FaceTopology;
    // typedef MonomialBasisProvider<dimension,Field> TestBasisProvider;
    // typedef MonomialBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    typedef OrthonormalBasisProvider<dimension,Field> TestBasisProvider;
    typedef OrthonormalBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef LagrangeBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LagrangeBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef LobattoBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LobattoBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef ScalarBasisProvider<dimension,Field> TestBasisProvider;
    // typedef ScalarBasisProvider<dimension-1,Field> TestFaceBasisProvider;

    typedef MonomialBasisProvider<dimension,Field> TestMBasisProvider;
    typedef MonomialBasisProvider<dimension-1,Field> TestFaceMBasisProvider;
    typedef SparseCoeffMatrix< Field, 1 > CoefficientMatrix;
    typedef StandardEvaluator< typename TestMBasisProvider::Basis > Evaluator;
    typedef PolynomialBasis< Evaluator, CoefficientMatrix > TestBasis;
    typedef StandardEvaluator< typename TestFaceMBasisProvider::Basis > FaceEvaluator;
    typedef PolynomialBasis< FaceEvaluator, CoefficientMatrix > TestFaceBasis;

    RaviartThomasL2Interpolation
      ( const unsigned int order,
      const FieldMatrix<Field,dimension+1,dimension> &normal )
      : order_(order),
        normal_(normal),
        mBasis_( TestBasisProvider::template basis<Topology>(order_-1) ),
        mFaceBasis_( TestFaceBasisProvider::template basis<FaceTopology>(order_) ),
        // mBasis_( TestBasisProvider::template basis<Topology>(TestBasisProvider::orthonormal,order_-1) ),
        // mFaceBasis_( TestFaceBasisProvider::template basis<FaceTopology>(TestFaceBasisProvider::orthonormal,order_) ),
        size_( (dimension+1)*mFaceBasis_.size()+dimension*mBasis_.size() )
    {}

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
        for (unsigned int i=0; i<mFaceBasis_.size(); ++i,++row)
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

      testBasisVal.resize(mFaceBasis_.size());
      testBasisInt.resize(mFaceBasis_.size());

      // mFaceBasis_.integral(testBasisInt);
      for (unsigned int f=0; f<dimension+1; ++f)
      {
        const typename SubQuadratureProvider::Quadrature &faceQuad = SubQuadratureProvider::template quadrature<Topology>( std::make_pair(f,2*order_+2) );
        const typename SubQuadratureProvider::SubQuadrature &faceSubQuad = SubQuadratureProvider::template subQuadrature<Topology>( std::make_pair(f,2*order_+2) );

        const unsigned int quadratureSize = faceQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          mFaceBasis_.template evaluate<0>(faceSubQuad.point(qi),testBasisVal);
          // for (int i=0;i<testBasisVal.size();++i)
          //  testBasisVal[i] /= testBasisInt[i];
          fillBnd( row, testBasisVal,
                   func.evaluate(faceQuad.point(qi)),
                   normal_[f], faceQuad.weight(qi),
                   func);
        }
        row += mFaceBasis_.size();
      }
      // element dofs
      if (row<size())
      {
        testBasisVal.resize(mBasis_.size());
        testBasisInt.resize(mBasis_.size());

        typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension, Field > QuadratureProvider;
        const typename QuadratureProvider::Quadrature &elemQuad = QuadratureProvider::template quadrature<Topology>(2*order_+1);
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
    FieldMatrix<Field,dimension+1,dimension> normal_;
    const TestBasis &mBasis_;
    const TestFaceBasis &mFaceBasis_;
    unsigned int size_;
  };

  template <class Topology, class RTInterpolation>
  struct RaviartThomasMatrix {
    typedef typename RTInterpolation::Field Field;
    enum {dim = Topology::dimension};
    struct VecMatrix
    {
      typedef MultiIndex<dim> MI;
      typedef MonomialBasis<Topology,MI> MIBasis;
      VecMatrix(unsigned int order)
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
      ~VecMatrix()
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
      // std::cout << matrix_ << std::endl;
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
    VecMatrix vecMatrix_;
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

    #define RTLAGRANGE 1
#if RTLAGRANGE
    typedef LagrangePointsCreator< ComputeField, dimension > PointsSetCreator;
    // typedef LobattoPointsCreator< ComputeField, dimension > PointsSetCreator;
    typedef RaviartThomasLagrangeInterpolation< ComputeField, PointsSetCreator > LocalInterpolation;
#else
    typedef RaviartThomasL2Interpolation< ComputeField, dimension > LocalInterpolation;
#endif

    template< class Topology >
    static const LocalInterpolation &localInterpolation ( const Key &order )
    {
      FieldMatrix<ComputeField,dimension+1,dimension> normal;
      for (int f=0; f<dimension+1; ++f)
        normal[f] = GenericGeometry::ReferenceElement<Topology,ComputeField>::integrationOuterNormal(f);
#if RTLAGRANGE
      const typename PointsSetCreator::LagrangePoints &points = PointsSetCreator::template lagrangePoints< Topology >( order+dimension );
      LocalInterpolation *interpolation = new LocalInterpolation(order,points,normal);
#else
      LocalInterpolation *interpolation = new LocalInterpolation(order,normal);
#endif
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
        /*
           const LocalCoefficients &keys = localCoefficients<Topology>(order);
           for (int index=0;index<keys.size();++index)
           std::cout << index << " -> "
                    << " (codim = " << keys.localKey(index).codim() << ", "
                    << "subentity = " << keys.localKey(index).subEntity() << ", "
                    << "index = " << keys.localKey(index).index() << "):" << std::endl;
         */
      }
      return *basis;
    }
    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      FieldMatrix<ComputeField,dimension+1,dimension> normal;
      for (int f=0; f<dimension+1; ++f)
        normal[f] = GenericGeometry::ReferenceElement<Topology,ComputeField>::integrationOuterNormal(f);
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
