// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASCOMPUTE_HH
#define DUNE_RAVIARTTHOMASCOMPUTE_HH
namespace Dune
{
  // A lagrange based interpolation for Raviart Thomas
  // --------------------------------------------------
  template< class F, class PointsSetCreator >
  class RaviartThomasLagrangeInterpolation
  {
    typedef typename PointsSetCreator::LagrangePoints LagrangePoints;
    static const unsigned int dimension = LagrangePoints::dimension;

    typedef RaviartThomasLagrangeInterpolation< F, PointsSetCreator > This;

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

    template< class Func >
    void interpolate ( Func &func ) const
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
  template< class F, unsigned int dimension >
  class RaviartThomasL2Interpolation
  {
    typedef RaviartThomasL2Interpolation< F, dimension > This;

  public:
    typedef F Field;
    typedef typename GenericGeometry::SimplexTopology<dimension>::type Topology;
    typedef typename GenericGeometry::SimplexTopology<dimension-1>::type FaceTopology;
    typedef MonomialBasisProvider<dimension,Field> TestBasisProvider;
    typedef MonomialBasisProvider<dimension-1,Field> TestFaceBasisProvider;

    RaviartThomasL2Interpolation
      ( const unsigned int order,
      const FieldMatrix<Field,dimension+1,dimension> &normal )
      : order_(order),
        normal_(normal),
        mBasis_( TestBasisProvider::template basis<Topology>(order_-1) ),
        mFaceBasis_( TestFaceBasisProvider::template basis<FaceTopology>(order_) ),
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


    void setLocalKeys(std::vector< LocalKey > &keys) const
    {
      keys.resize(size());
      unsigned int row = 0;
      for (int f=0; f<dimension+1; ++f)
        for (int i=0; i<mFaceBasis_.size(); ++i,++row)
          keys[row] = LocalKey(f,1,i);
      for (int i=0; i<mBasis_.size()*dimension; ++i,++row)
        keys[row] = LocalKey(0,0,i);
      assert( row == size() );
    }

    template< class Func >
    void interpolate ( Func &func ) const
    {
      std::vector< Field > testBasisVal;

      for (unsigned int i=0; i<size(); ++i)
        for (unsigned int j=0; j<func.size(); ++j)
          func.set(i,j,0);

      unsigned int row = 0;

      // boundary dofs:
      typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension-1, Field > FaceQuadratureProvider;
      typedef Dune::GenericGeometry::SubQuadratureProvider< dimension, FaceQuadratureProvider> SubQuadratureProvider;

      testBasisVal.resize(mFaceBasis_.size());

      for (unsigned int f=0; f<dimension+1; ++f)
      {
        const typename SubQuadratureProvider::Quadrature &faceQuad = SubQuadratureProvider::template quadrature<Topology>( std::make_pair(f,2*order_+1) );
        const typename SubQuadratureProvider::SubQuadrature &faceSubQuad = SubQuadratureProvider::template subQuadrature<Topology>( std::make_pair(f,2*order_+1) );

        const unsigned int quadratureSize = faceQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          mFaceBasis_.template evaluate<0>(faceSubQuad.point(qi),testBasisVal);
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

        typedef Dune::GenericGeometry::GenericQuadratureProvider< dimension, Field > QuadratureProvider;
        const typename QuadratureProvider::Quadrature &elemQuad = QuadratureProvider::template quadrature<Topology>(2*order_+1);
        const unsigned int quadratureSize = elemQuad.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          mBasis_.template evaluate<0>(elemQuad.point(qi),testBasisVal);
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
          matrix.add(row,col, weight*(cFactor*(*miter)) );
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
            matrix.add(row+i,col, weight*(*rtiter)[i]*(*miter) );
          }
        }
        assert( miter == mVal.end() );
      }
    }

    unsigned int order_;
    FieldMatrix<Field,dimension+1,dimension> normal_;
    const typename TestBasisProvider::Basis &mBasis_;
    const typename TestFaceBasisProvider::Basis &mFaceBasis_;
    unsigned int size_;
  };

  template <class Topology, class Field>
  struct RTPreMatrix
  {
    enum {dim = Topology::dimension};
    typedef MultiIndex<dim> MI;
    typedef MonomialBasis<Topology,MI> MIBasis;
    RTPreMatrix(unsigned int order)
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
          unsigned int q = -1;
          MI xval = val[notHomogen+i];
          xval *= x[r];
          for (unsigned int w=homogen+notHomogen; w<val.size(); ++w)
          {
            if (val[w] == xval)
            {
              q=w;
              break;
            }
          }
          assert(q>=0);
          mat_[row][q] = 1.;
          ++row;
        }
      }
    }
    ~RTPreMatrix()
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

}
#endif // DUNE_RAVIARTTHOMASBASIS_HH
