// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASPREBASIS_HH
#define DUNE_RAVIARTTHOMASPREBASIS_HH

#include <fstream>
#include <utility>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/polynomialbasis.hh>

namespace Dune
{
  template < class Topology, class Field >
  struct RTVecMatrix;

  template <unsigned int dim, class Field>
  struct RTPreBasisFactory
  {
    typedef MonomialBasisProvider<dim,Field> MBasisFactory;
    typedef typename MBasisFactory::Object MBasis;
    typedef StandardEvaluator<MBasis> EvalMBasis;
    typedef PolynomialBasisWithMatrix<EvalMBasis,SparseCoeffMatrix<Field,dim> > Basis;

    typedef const Basis Object;
    typedef std::size_t Key;

    template <unsigned int dd, class FF>
    struct EvaluationBasisFactory
    {
      typedef MonomialBasisProvider<dd,FF> Type;
    };
    template< class Topology >
    static Object *create ( const Key &order )
    {
      RTVecMatrix<Topology,Field> vecMatrix(order);
      MBasis *mbasis = MBasisFactory::template create<Topology>(order+1);
      typename std::remove_const<Object>::type *tmBasis =
        new typename std::remove_const<Object>::type(*mbasis);
      tmBasis->fill(vecMatrix);
      return tmBasis;
    }
    static void release( Object *object ) { delete object; }
  };
  template <class Topology, class Field>
  struct RTVecMatrix
  {
    static const unsigned int dim = Topology::dimension;
    typedef MultiIndex<dim,Field> MI;
    typedef MonomialBasis<Topology,MI> MIBasis;
    RTVecMatrix(std::size_t order)
    {
      MIBasis basis(order+1);
      FieldVector< MI, dim > x;
      for( unsigned int i = 0; i < dim; ++i )
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
      for (unsigned int i=0; i<rows(); ++i) {
        delete [] mat_[i];
      }
      delete [] mat_;
    }
    unsigned int cols() const {
      return col_;
    }
    unsigned int rows() const {
      return row_;
    }
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      const unsigned int N = cols();
      assert( vec.size() == N );
      for (unsigned int i=0; i<N; ++i)
        field_cast(mat_[row][i],vec[i]);
    }
    unsigned int row_,col_,row1_;
    Field **mat_;
  };


}
#endif // DUNE_RAVIARTTHOMASPREBASIS_HH
