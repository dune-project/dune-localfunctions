// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASBASIS_HH
#define DUNE_RAVIARTTHOMASBASIS_HH
#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>

#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/polynomialbasis.hh>
namespace Dune
{
  template< class Topology, class F >
  class RaviartThomasLagrangeInterpolation
  {
    typedef RaviartThomasLagrangeInterpolation< Topology, F > This;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

  private:
    typedef Dune::LagrangePoints< Field, dimension > LagrangePoints;
    typedef Dune::LagrangePointsCreator< Field, dimension > LagrangePointsCreator;

  public:
    RaviartThomasLagrangeInterpolation ( const unsigned int order )
      : lagrangePoints_ (LagrangePointsCreator::template lagrangePoints< Topology >( order+dimension ))
    {}

    template< class Eval, class Matrix >
    void interpolate ( Eval &eval, Matrix &coefficients )
    {
      coefficients.resize( eval.size(), eval.size( ) );
      typedef typename LagrangePoints::iterator Iterator;
      unsigned int row = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
      {
        if (it->localKey().codim()>1)
          continue;
        else if (it->localKey().codim()==1) {
          fillBnd ( row, it->localKey(), eval.template evaluate<0>(it->point()),
                    coefficients );
        }
        else
          fillInterior ( row, it->localKey(), eval.template evaluate<0>(it->point()),
                         coefficients );
      }
    }

  private:
    /** /brief evaluate boundary functionals **/
    template <class LocalKey, class Iterator,class Matrix>
    void fillBnd (unsigned int &row,
                  const LocalKey &key, Iterator iter, Matrix &matrix) const
    {
      const Dune::FieldVector<double,dimension> &normal = GenericGeometry::ReferenceElement<Topology,double>::integrationOuterNormal(key.subEntity());
      unsigned int col = 0;
      for ( ; !iter.done() ; ++iter,++col) {
        matrix(row,col) = 0.;
        for (unsigned int d=0; d<dimension; ++d) {
          matrix(row,col) += iter->block()[d]*normal[d];
        }
      }
      ++row;
    }
    template <class LocalKey, class Iterator,class Matrix>
    void fillInterior (unsigned int &row,
                       const LocalKey &key, Iterator iter,Matrix &matrix) const
    {
      unsigned int col = 0;
      for ( ; !iter.done() ; ++iter,++col)
        for (unsigned int d=0; d<dimension; ++d)
          matrix(row+d,col) = iter->block()[d];
      row+=dimension;
    }
    const LagrangePoints &lagrangePoints_;
  };

#if 0
  template< class Topology, class F >
  class RaviartThomasL2Interpolation
  {
    typedef RaviartThomasL2Interpolation< Topology, F > This;

  public:
    typedef F Field;
    static const unsigned int dimension = Topology::dimension;

    RaviartThomasLagrangeInterpolation ( const unsigned int order )
    {}

    template< class Eval, class Matrix >
    void interpolate ( Eval &eval, Matrix &coefficients )
    {
      coefficients.resize( eval.size(), eval.size( ) );
      for (int i=0; i<eval.size(); ++i)
        for (int j=0; j<eval.size(); ++j)
          coefficients(i,j) = Zero<Field>();
      unsigned int row;

      unsigned int nrFaces = GenericGeometry::ReferenceElement<Topology>::template Codim<1>::size;

      // boundary dofs:
      for (int f=0; f<nrFaces; ++f)
      {
        // get subentity info (topology for quadrature and basis)
        // and the hybrid mapping to the codim<0> entity
        const Dune::FieldVector<double,dimension> &normal = GenericGeometry::ReferenceElement<Topology,double>::integrationOuterNormal(key.subEntity());
        GenericQuadrature< dimension-1, Field >
        quadrature( subId, 2*order_+1 );
        typedef MonomialBasisProvider<dimension-1,Field> MBasisProvider;
        const typename MBasisProvider::Basis &mBasis = BasisProvider::basis(subId,order_);
        StandardEvaluator<typename MBasisProvider::Basis> mEval(mBasis);
        const unsigned int quadratureSize = quadrature.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          const Dune::FieldVector<double,dimension> point =
            faceMap.global(quadrature.point(qi));
          fillBnd( row,
                   mEval.template evaluate<0>(quadrature.point(qi)),
                   eval.template evaluate<0>(point),
                   normal,
                   quadrature.weight(qi),
                   coefficients);
        }
        row += mEval.size();
      }
      // element dofs
      if (row<eval.size())
      {
        GenericQuadrature< dimension, Field >
        quadrature( id, 2*order_+1 );
        typedef MonomialBasisProvider<dimension,Field> MBasisProvider;
        const typename MBasisProvider::Basis &mBasis = BasisProvider::basis(id,order_-1);
        StandardEvaluator<typename MBasisProvider::Basis> mEval(mBasis);
        const unsigned int quadratureSize = quadrature.size();
        for( unsigned int qi = 0; qi < quadratureSize; ++qi )
        {
          fillInterior( row,
                        mEval.template evaluate<0>(quadrature.point(qi)),
                        eval.template evaluate<0>(quadrature.point(qi)),
                        quadrature.weight(qi),
                        coefficients );
        }
        row += mEval.size()*dimension;
      }
      assert(row==eval.size());
    }

  private:
    /** /brief evaluate boundary functionals **/
    template <class MIterator, class RTIterator,class Matrix>
    void fillBnd (unsigned int startRow,
                  rowIterator riter,
                  colIterator citer,
                  const FieldVector<double,dimension> &normal,
                  double weight,
                  Matrix &matrix) const
    {
      unsigned int col = 0;
      for ( ; !citer.done() ; ++citer,++col)
      {
        Field cFactor = Zero<Field>();
        for (unsigned int d=0; d<dimension; ++d) {
          cFactor += citer->block()[d]*normal[d];
          for (unsigned int row = startRow ; !riter.done(); ++riter; ++row )
          {
            matrix(row,col) += weight*(cFactor*riter->block()[0]);
          }
        }
      }
      template <class MIterator, class RTIterator,class Matrix>
      void fillBnd (unsigned int startRow,
                    rowIterator riter,
                    colIterator citer,
                    double weight,
                    Matrix &matrix) const
      {
        unsigned int col = 0;
        for ( ; !citer.done() ; ++citer,++col)
        {
          for (unsigned int row = startRow ; !riter.done(); ++riter; )
          {
            for (int i=0; i<dimension; ++i,++row)
            {
              matrix(row,col) += weight*(citer->block()[i]*riter->block()[0]);
            }
          }
        }
      }
    };
#endif

  template <class Topology,class scalar_t>
  struct RaviartThomasMatrix {
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
        mat_ = new double*[row_];
        int row = 0;
        for (unsigned int i=0; i<notHomogen+homogen; ++i)
        {
          for (unsigned int r=0; r<dim; ++r)
          {
            for (unsigned int rr=0; rr<dim; ++rr)
            {
              mat_[row] = new double[col_];
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
            mat_[row] = new double[col_];
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
      const scalar_t operator() ( int r, int c ) const
      {
        assert(r<(int)row_ && c<(int)col_);
        return mat_[r][c];
      }
      unsigned int row_,col_,row1_;
      double **mat_;
    };
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    typedef MonomialBasis<Topology,scalar_t> MBasis;
    RaviartThomasMatrix(unsigned int order) : order_(order),
                                              vecMatrix_(order)
    {
      MBasis basis(order+1);
      typedef StandardEvaluator<MBasis> EvalMBasis;
      typedef PolynomialBasisWithMatrix<EvalMBasis,SparseCoeffMatrix<scalar_t,dim> > TMBasis;
      TMBasis tmBasis(basis);
      tmBasis.fill(vecMatrix_);
      StandardEvaluator<TMBasis> tmEval(tmBasis);
      RaviartThomasLagrangeInterpolation< Topology, scalar_t  > interpolation(order_);
      interpolation.interpolate( tmEval , matrix_ );
      matrix_.invert();
    }
    unsigned int colSize(int row) const {
      return vecMatrix_.colSize(row);
    }
    unsigned int rowSize() const {
      return vecMatrix_.rowSize();
    }
    const scalar_t operator() ( int r, int c ) const
    {
      int rmod = r%dim;
      int rdiv = r/dim;
      scalar_t ret = 0;
      for (unsigned int k=0; k<matrix_.cols(); ++k) {
        ret += matrix_(k,rdiv)*vecMatrix_(k*dim+rmod,c);
      }
      return ret;
    }
    int order_;
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
    typedef unsigned int Key;
    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    template <class Topology>
    static Basis &basis(unsigned int order)
    {
      const MBasis &_mBasis = MonomialBasisProvider<dimension,StorageField>::template basis<SimplexTopology>(order+1);
      Basis *basis = new Basis(_mBasis);
      RaviartThomasMatrix<Topology,ComputeField> matrix(order);
      basis->fill(matrix);
      /*
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
       */
      return *basis;
    }
    static void release ( const Basis &basis )
    {
      delete &basis;
    }
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct RaviartThomasBasisProvider
    : public BasisProvider<RaviartThomasBasisCreator<dim,SF,CF> >
  {};
}
#endif // DUNE_RAVIARTTHOMASBASIS_HH
