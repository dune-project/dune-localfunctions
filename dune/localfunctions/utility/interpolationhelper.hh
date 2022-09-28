// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef GENERIC_INTERPOLATIONHELPER_HH
#define GENERIC_INTERPOLATIONHELPER_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/concept.hh>
#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  // A small helper class to avoid having to
  // write the interpolation twice (once for function
  // and once for a basis)
  template< class F, unsigned int dimension >
  struct InterpolationHelper
  {
    template <class Func,class Container, bool type>
    struct Helper;
  };
  template <class F,unsigned int d>
  template <class Func,class Vector>
  struct InterpolationHelper<F,d>::Helper<Func,Vector,true>
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
    template <class DomainVector,
              std::enable_if_t<models<Impl::FunctionWithCallOperator<DomainVector>, Func>(), int> = 0>
    const Result &evaluate(const DomainVector &x) const
    {
      field_cast(func_(x), tmp_[0] );
      return tmp_;
    }
    template <class DomainVector,
              std::enable_if_t<not models<Impl::FunctionWithCallOperator<DomainVector>, Func>(), int> = 0>
    const Result &evaluate(const DomainVector &x) const
    {
      typename Func::DomainType xx ;
      typename Func::RangeType ff ;
      field_cast(x,xx);
      func_.evaluate(xx,ff);
      field_cast(ff, tmp_[0] );
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
  struct InterpolationHelper<F,d>::Helper<Basis,Matrix,false>
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
}
#endif // GENERIC_INTERPOLATIONHELPER_HH
