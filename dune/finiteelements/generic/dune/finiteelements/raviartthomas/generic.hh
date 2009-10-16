// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_HH
#define DUNE_GENERIC_HH
#include <fstream>
#include <utility>
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

  template< unsigned int dimension, class Impl >
  struct InterpolationWrapper
  {
    // A small helper class to avoid having to
    // write the interpolation twice (once for function
    // and once for a basis)
    template <class F, class Func,class Container, bool type>
    struct Helper;

    InterpolationWrapper(const Impl *impl)
      : impl_(impl) {}
    ~InterpolationWrapper()
    {
      delete impl_;
    }

    unsigned int order() const
    {
      return impl_->order();
    }
    unsigned int size() const
    {
      return impl_->size();
    }

    template< class Function, class Fy >
    void interpolate ( const Function &function, std::vector< Fy > &coefficients ) const
    {
      coefficients.resize(size());
      Helper<typename Impl::Field,Function,std::vector<Fy>,true> func( function, coefficients );
      impl_->interpolate(func);
    }
    template< class Basis, class Matrix >
    void interpolate ( const Basis &basis, Matrix &matrix ) const
    {
      matrix.resize( size(), basis.size() );
      Helper<typename Impl::Field,Basis,Matrix,false> func( basis, matrix );
      impl_->interpolate(func);
    }

  private:
    const Impl *impl_;
  };
  template <unsigned int d,class Impl>
  template <class F, class Func,class Vector>
  struct InterpolationWrapper<d,Impl>::Helper<F,Func,Vector,true>
  // Func is of Function type
  {
    typedef Dune::FieldVector<F,d> Result[1];
    Helper(const Func & func, Vector &vec)
      : func_(func),
        vec_(vec)
    {}
    const typename Vector::value_type &operator()(unsigned int row,unsigned int col)
    {
      return vec_[row];
    }
    template <class Fy>
    void set(unsigned int row,unsigned int col,
             const Fy &val)
    {
      assert(col==0);
      assert(row<vec_.size());
      field_cast( val, vec_[row] );
    }
    template <class Fy>
    void add(unsigned int row,unsigned int col,
             const Fy &val)
    {
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
  template <unsigned int d,class Impl>
  template <class F,class Basis,class Matrix>
  struct InterpolationWrapper<d,Impl>::Helper<F,Basis,Matrix,false>
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

  // ********************************************************

  template <class Topology,
      class PreMatrix>
  struct MatrixBuilder {
    enum {dim = Topology::dimension};
    typedef typename PreMatrix::Field Field;
    typedef Dune::AlgLib::Matrix< Field > mat_t;
    template <class Interpolation,class Basis>
    MatrixBuilder(const Interpolation &interpolation,
                  const Basis &basis,
                  const PreMatrix &preMatrix) :
      preMatrix_(preMatrix)
    {
      typedef StandardEvaluator<Basis> EvalMBasis;
      typedef PolynomialBasisWithMatrix<EvalMBasis,SparseCoeffMatrix<Field,dim> > TMBasis;
      TMBasis tmBasis(basis);
      tmBasis.fill(preMatrix_);
      interpolation.interpolate( tmBasis , matrix_ );
      matrix_.invert();
    }
    unsigned int colSize(int row) const {
      return preMatrix_.colSize(row);
    }
    unsigned int rowSize() const {
      return preMatrix_.rowSize();
    }
    const Field operator() ( int r, int c ) const
    {
      int rmod = r%dim;
      int rdiv = r/dim;
      Field ret = 0;
      for (unsigned int k=0; k<matrix_.cols(); ++k) {
        ret += matrix_(k,rdiv)*preMatrix_(k*dim+rmod,c);
      }
      return ret;
    }
    const PreMatrix &preMatrix_;
    mat_t matrix_;
  };


}
#endif // DUNE_RAVIARTTHOMASBASIS_HH
