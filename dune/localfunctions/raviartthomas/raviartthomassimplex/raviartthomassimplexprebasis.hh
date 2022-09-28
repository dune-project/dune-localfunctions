// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RAVIARTTHOMASPREBASIS_HH
#define DUNE_RAVIARTTHOMASPREBASIS_HH

#include <fstream>
#include <utility>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/polynomialbasis.hh>

namespace Dune
{
  template < GeometryType::Id geometryId, class Field >
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
    template< GeometryType::Id geometryId >
    static Object *create ( const Key &order )
    {
      RTVecMatrix<geometryId,Field> vecMatrix(order);
      MBasis *mbasis = MBasisFactory::template create<geometryId>(order+1);
      typename std::remove_const<Object>::type *tmBasis = new typename std::remove_const<Object>::type(*mbasis);
      tmBasis->fill(vecMatrix);
      return tmBasis;
    }
    static void release( Object *object ) { delete object; }
  };

  template <GeometryType::Id geometryId, class Field>
  struct RTVecMatrix
  {
    static constexpr GeometryType geometry = geometryId;
    static const unsigned int dim = geometry.dim();
    typedef MultiIndex<dim,Field> MI;
    typedef MonomialBasis<geometryId,MI> MIBasis;
    RTVecMatrix(std::size_t order)
    {
      /*
       * Construction of Raviart-Thomas elements in high dimensions see "Mixed Finite Elements in \R^3" by Nedelec, 1980.
       *
       * Let $\P_{n,k}$ be the space of polynomials in $n$ variables with degree $\leq k$.
       * The space of Raviart-Thomas functions in $n$ dimensions with index $k$ is defined as
       *
       * \begin{equation*}
       *    RT_k := (\P_{k-1})^n \oplus \widetilde \P_k x
       * \end{equation*}
       * with $x=(x_1,x_2,\dots, x_n)$ in $n$ dimensions and $\widetilde \P_k$ the homogeneous polynomials of degree $k$.
       *
       * For $RT_k$ holds
       * \begin{equation*}
       *    (\P_{k-1})^n \subset RT_k \subset (\P_k)^n.
       * \end{equation*}
       *
       * We construct $(\P_k)^n$ and and only use the monomials contained in $RT_k$.
       *
       */

      MIBasis basis(order+1);
      FieldVector< MI, dim > x;
      /*
       * Init MultiIndices
       * x[0]=(1,0,0) x
       * x[1]=(0,1,0) y
       * x[2]=(0,0,1) z
       */
      for( unsigned int i = 0; i < dim; ++i )
        x[ i ].set( i, 1 );
      std::vector< MI > val( basis.size() );

      // val now contains all monomials in $n$ dimensions with degree $\leq order+1$
      basis.evaluate( x, val );

      col_ = basis.size();

      // get $\dim (\P_{order-1})$
      unsigned int notHomogen = 0;
      if (order>0)
        notHomogen = basis.sizes()[order-1];

      // get $\dim \widetilde (\P_order)$
      unsigned int homogen = basis.sizes()[order]-notHomogen;

      /*
       *
       * The set $RT_k$ is defined as
       *
       * \begin{equation}
       *   RT_k :=  (\P_k)^dim  +  \widetilde \P_k x   with x\in \R^n.
       * \end{equation}
       *
       * The space $\P_k$ is split in $\P_k = \P_{k-1} + \widetilde \P_k$.
       *
       * \begin{align}
       *   RT_k &=  (\P_{k-1} + \widetilde \P_k)^dim  +  \widetilde \P_k x   with x\in \R^n
       *        &=  (\P_{k-1})^n + (\widetilde \P_k)^n  +  \widetilde \P_k x   with x\in \R^n
       * \end{align}
       *
       * Thus $\dim RT_k = n * \dim \P_{k-1} + (n+1)*\dim (\widetilde \P_k)$
       */

      // row_ = \dim RT_k *dim
      row_ = (notHomogen*dim+homogen*(dim+1))*dim;
      mat_ = new Field*[row_];
      int row = 0;

      /* Assign the correct values for the nonhomogeneous polymonials $p\in (\P_{oder-1})^dim$
       * A basis function is represented by $dim$ rows.
       */
      for (unsigned int i=0; i<notHomogen+homogen; ++i)
      {
        for (unsigned int r=0; r<dim; ++r)
        {
          for (unsigned int rr=0; rr<dim; ++rr)
          {
            // init row to zero
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

      /* Assign the correct values for the homogeneous polymonials $p\in RT_k \backslash (\P_{oder-1})^dim$
       * A basis function is represented by $dim$ rows.
       */
      for (unsigned int i=0; i<homogen; ++i)
      {
        for (unsigned int r=0; r<dim; ++r)
        {
          // init rows to zero
          mat_[row] = new Field[col_];
          for (unsigned int j=0; j<col_; ++j)
          {
            mat_[row][j] = 0.;
          }

          unsigned int w;
          // get a monomial $ p \in \widetilde \P_{order}$
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
    unsigned int row_,col_;
    Field **mat_;
  };


}
#endif // DUNE_RAVIARTTHOMASPREBASIS_HH
