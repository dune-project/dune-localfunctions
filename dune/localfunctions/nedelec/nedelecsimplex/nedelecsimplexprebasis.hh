// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXPREBASIS_HH
#define DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXPREBASIS_HH

#include <fstream>
#include <utility>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/polynomialbasis.hh>

namespace Dune
{
  template < GeometryType::Id geometryId, class Field >
  struct NedelecVecMatrix;

  template <unsigned int dim, class Field>
  struct NedelecPreBasisFactory
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
    static Object *create ( Key order )
    {
      /*
       * The nedelec parameter begins at 1.
       * This is the numbering used by J.C. Nedelec himself.
       * See "Mixed Finite Elements in \R^3" published in 1980.
       *
       * This construction is based on the construction of Raviart-Thomas elements.
       * There the numbering starts at 0.
       * Because of this we reduce the order internally by 1.
       */
      order--;
      NedelecVecMatrix<geometryId,Field> vecMatrix(order);
      MBasis *mbasis = MBasisFactory::template create<geometryId>(order+1);
      std::remove_const_t<Object>* tmBasis = new std::remove_const_t<Object>(*mbasis);
      tmBasis->fill(vecMatrix);
      return tmBasis;
    }
    static void release( Object *object ) { delete object; }
  };

  template <GeometryType::Id geometryId, class Field>
  struct NedelecVecMatrix
  {
    static constexpr GeometryType geometry = geometryId;
    static const unsigned int dim = geometry.dim();
    typedef MultiIndex<dim,Field> MI;
    typedef MonomialBasis<geometryId,MI> MIBasis;
    NedelecVecMatrix(std::size_t order)
    {
      /*
       * Construction of Nedelec elements see "Mixed Finite Elements in \R^3" by Nedelec, 1980.
       *
       * Let $\P_{n,k}$ be the space of polynomials in $n$ variables with degree $\leq k$.
       * The space of Nedelec functions in $n$ dimensions with index $k$ is defined as
       *
       * \begin{equation*}
       *    Ned_k := (\P_{n,k-1})^n \oplus \{p \in (\P_{n,k})^n: <p,x>=0 \}
       * \end{equation*}
       * with $x=(x,y)$ in two dimensions and $x=(x,y,z)$ in three dimensions.
       *
       * For $Ned_k$ holds
       * \begin{equation*}
       *    (\P_{n,k-1})^n \subset Ned_k \subset (\P_{n,k})^n.
       * \end{equation*}
       *
       * We construct $(\P_{n,k})^n$ and and only use the monomials contained in $Ned$.
       *
       */
      if( (dim!=2 && dim!=3) || !geometry.isSimplex())
        DUNE_THROW(Dune::NotImplemented,"High order nedelec elements are only supported by simplices in 2d and 3d");

      MIBasis basis(order+1);
      FieldVector< MI, dim > x;
      /*
       * Init MultiIndices
       * x[0]=(1,0,0) x
       * x[1]=(0,1,0) y
       * x[2]=(0,0,1) z
       */
      for( unsigned int i = 0; i < dim; ++i )
        x[i].set(i,1);
      std::vector< MI > val( basis.size() );

      // val now contains all monomials in $n$ dimensions with degree $\leq order+1$
      basis.evaluate( x, val );

      col_ = basis.size();

      // get $\dim (\P_{n,order-1})$
      unsigned int notHomogen = 0;
      if (order>0)
        notHomogen = basis.sizes()[order-1];

      // the number of basis functions for the set of homogeneous polynomials with degree $order$
      unsigned int homogen = basis.sizes()[order]-notHomogen;

      /*
       * 2D:
       * \begin{equation*}
       *    Ned_{order} =  (\P_{order-1})^2 \oplus (-y,x)^T \widetilde \P_{order-1}
       * \end{equation*}
       *
       * It gets more complicated in higher dimensions.
       *
       * 3D:
       * \begin{equation*}
       *    Ned_{order} =  (\P_{n,order-1})^3 \oplus (z,0,-x)^T \widetilde \P_{n,order-1} \oplus (-y,x,0)^T \widetilde \P_{n,order-1} \oplus (0,-z,y)^T \widetilde \P_{n-1,order-1}
       * \end{equation*}
       *
       * Note the last term. The index $n-1$ is on purpose.
       * Else i.e. k=2
       *
       *  (0,z,-y)^T x = (z,0,-x)^T y - (y,-x,0)^T z.
       *
       */

      /*
       * compute the number of rows for the coefficient matrix
       *
       * row_ = dim* \dim Ned_{order}
       */
      if (dim == 2)
        row_ = (notHomogen*dim+homogen*(dim+1))*dim;
      else if (dim==3)
      {
        // get dim \P_{n-1,order-1}
        int homogenTwoVariables = 0;
        for( int w = notHomogen; w<notHomogen + homogen; w++)
          if (val[w].z(0)==0)
            homogenTwoVariables++;
        row_ = (notHomogen*dim+homogen*(dim+2) + homogenTwoVariables)*dim;
      }

      mat_ = new Field*[row_];
      int row = 0;

      /* Assign the correct values for the nonhomogeneous polymonials $p\in (\P_{n,order-1})^dim$
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
              mat_[row][j] = 0.;

            if (r==rr)
              mat_[row][i] = 1.;
            ++row;
          }
        }
      }

      /* Assign the correct values for the homogeneous polymonials $p\in Ned_{order} \backslash (\P_{n,order-1})^dim$
       * A basis function is represented by $dim$ rows.
       */
      for (unsigned int i=0; i<homogen; ++i)
      {
        // get a monomial $ p \in \P_{n,order}\backslash \P_{n,order-1}$
        MI xval = val[notHomogen+i];
        if(dim==2)
        {
          for (unsigned int r=0; r<dim; ++r)
          {
            // init rows to zero
            mat_[row+r] = new Field[col_];
            for (unsigned int j=0; j<col_; ++j)
              mat_[row+r][j] = 0.;
          }

          /* set $(-y,x)^T p$ with a homogeneous monomial $p$
           *
           * The loop over the monomials is needed to obtain the corresponding column index.
           */
          for (int w=homogen+notHomogen; w<val.size(); ++w)
          {
            if (val[w] == xval*x[0])
              mat_[row+1][w] = 1.;
            if (val[w] == xval*x[1])
              mat_[row][w] = -1.;
          }
          row +=dim;
        }
        else if(dim==3)
        {
          int skipLastDim = xval.z(0)>0;
          /*
           * Init 9 rows to zero.
           * If the homogeneous monomial has a positive x-exponent (0,-z,y) gets skipped ( see example for the Nedelec space in 3D )
           * In this case only 6 rows get initialised.
           */
          for (unsigned int r=0; r<dim*(dim-skipLastDim); ++r)
          {
            // init rows to zero
            mat_[row+r] = new Field[col_];
            for (unsigned int j=0; j<col_; ++j)
              mat_[row+r][j] = 0.;
          }

          /*
           * first $dim$ rows are for (z,0,-x)
           *
           * second $dim$ rows are for (-y,x,0)
           *
           * third $dim$ rows are for (0,-z,y)
           *
           */
          for (unsigned int r=0; r<dim - skipLastDim; ++r)
          {
            int index = (r+dim-1)%dim;
            for (int w=homogen+notHomogen; w<val.size(); ++w)
            {
              if (val[w] == xval*x[index])
                mat_[row+r][w] = 1.;
              if (val[w] == xval*x[r])
                mat_[row+index][w] = -1.;
            }
            row +=dim;
          }

        }
      }
    }

    ~NedelecVecMatrix()
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
#endif // #ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXPREBASIS_HH
