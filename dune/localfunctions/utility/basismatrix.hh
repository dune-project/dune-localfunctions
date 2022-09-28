// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_BASISMATRIX_HH
#define DUNE_BASISMATRIX_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/utility/lfematrix.hh>
#include <dune/localfunctions/utility/monomialbasis.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>

namespace Dune
{
  /****************************************
  * A dense matrix representation of a ''polynomial''
  * basis. Its represent a basis as a linear
  * combination of a second basis, i.e., a
  * monomial basis. It is simular to the PolynomialBasis
  * but it not derived from the LocalBasis class.
  * It is used to define a ''pre basis''.
  ****************************************/
  template< class PreBasis, class Interpolation,
      class Field >
  struct BasisMatrix;

  template< class PreBasis, class Interpolation,
      class Field >
  struct BasisMatrixBase : public LFEMatrix<Field>
  {
    typedef LFEMatrix<Field> Matrix;

    BasisMatrixBase( const PreBasis& preBasis,
                     const Interpolation& localInterpolation )
      : cols_(preBasis.size())
    {
      localInterpolation.interpolate( preBasis, *this );

      if ( !Matrix::invert() )
      {
        DUNE_THROW(MathError, "While computing basis a singular matrix was constructed!");
      }
    }
    unsigned int cols () const
    {
      return cols_;
    }
    unsigned int rows () const
    {
      return Matrix::rows();
    }
  private:
    unsigned int cols_;
  };

  template< GeometryType::Id geometryId, class F,
      class Interpolation,
      class Field >
  struct BasisMatrix< const MonomialBasis< geometryId, F >, Interpolation, Field >
    : public BasisMatrixBase< const MonomialBasis< geometryId, F >, Interpolation, Field >
  {
    typedef const MonomialBasis< geometryId, F > PreBasis;
    typedef BasisMatrixBase<PreBasis,Interpolation,Field> Base;
    typedef typename Base::Matrix Matrix;

    BasisMatrix( const PreBasis& preBasis,
                 const Interpolation& localInterpolation )
      : Base(preBasis, localInterpolation)
    {}
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      const unsigned int N = Matrix::rows();
      assert( Matrix::cols() == N && vec.size() == N );
      // note: that the transposed matrix is computed,
      //       and is square
      for (unsigned int i=0; i<N; ++i)
        field_cast(Matrix::operator()(i,row),vec[i]);
    }
  };
  template< int dim, class F,
      class Interpolation,
      class Field >
  struct BasisMatrix< const Dune::VirtualMonomialBasis< dim, F >, Interpolation, Field >
    : public BasisMatrixBase< const VirtualMonomialBasis< dim, F >, Interpolation, Field >
  {
    typedef const VirtualMonomialBasis< dim, F > PreBasis;
    typedef BasisMatrixBase<PreBasis,Interpolation,Field> Base;
    typedef typename Base::Matrix Matrix;

    BasisMatrix( const PreBasis& preBasis,
                 const Interpolation& localInterpolation )
      : Base(preBasis, localInterpolation)
    {}
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      const unsigned int N = Matrix::rows();
      assert( Matrix::cols() == N && vec.size() == N );
      // note: that the transposed matrix is computed,
      //       and is square
      for (unsigned int i=0; i<N; ++i)
        field_cast(Matrix::operator()(i,row),vec[i]);
    }
  };
  template< class Eval, class CM, class D, class R,
      class Interpolation,
      class Field >
  struct BasisMatrix< const PolynomialBasis<Eval,CM,D,R>, Interpolation, Field >
    : public BasisMatrixBase< const PolynomialBasis<Eval,CM,D,R>, Interpolation, Field >
  {
    typedef const PolynomialBasis<Eval,CM,D,R> PreBasis;
    typedef BasisMatrixBase<PreBasis,Interpolation,Field> Base;
    typedef typename Base::Matrix Matrix;

    BasisMatrix( const PreBasis& preBasis,
                 const Interpolation& localInterpolation )
      : Base(preBasis, localInterpolation),
        preBasis_(preBasis)
    {}
    unsigned int cols() const
    {
      return preBasis_.matrix().baseSize() ;
    }
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      assert( Matrix::rows() == Matrix::cols() );
      assert( vec.size() == preBasis_.matrix().baseSize() );
      assert( Matrix::cols() == preBasis_.size() );
      for (unsigned int j=0; j<Matrix::cols(); ++j)
        vec[j] = 0;
      for (unsigned int i=0; i<Matrix::rows(); ++i)
        preBasis_.matrix().
        addRow(i,Base::Matrix::operator()(i,row),vec);
    }
  private:
    const PreBasis& preBasis_;
  };
  template< class Eval, class CM,
      class Interpolation,
      class Field >
  struct BasisMatrix< const PolynomialBasisWithMatrix<Eval,CM>, Interpolation, Field >
    : public BasisMatrixBase< const PolynomialBasisWithMatrix<Eval,CM>, Interpolation, Field >
  {
    typedef const PolynomialBasisWithMatrix<Eval,CM> PreBasis;
    typedef BasisMatrixBase<PreBasis,Interpolation,Field> Base;
    typedef typename Base::Matrix Matrix;

    BasisMatrix( const PreBasis& preBasis,
                 const Interpolation& localInterpolation )
      : Base(preBasis, localInterpolation),
        preBasis_(preBasis)
    {}
    unsigned int cols() const
    {
      return preBasis_.matrix().baseSize() ;
    }
    unsigned int rows () const
    {
      assert( Matrix::rows() == preBasis_.matrix().size() );
      return preBasis_.matrix().size()*CM::blockSize ;
    }
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      unsigned int r = row / CM::blockSize;
      assert( r < Matrix::rows() );
      assert( Matrix::rows() == Matrix::cols() );
      assert( vec.size() == preBasis_.matrix().baseSize() );
      assert( Matrix::cols() == preBasis_.size() );
      for (unsigned int j=0; j<vec.size(); ++j)
        vec[j] = 0;
      for (unsigned int i=0; i<Matrix::rows(); ++i)
        preBasis_.matrix().
        addRow(i*CM::blockSize+row%CM::blockSize,Base::Matrix::operator()(i,r),vec);
    }
  private:
    const PreBasis& preBasis_;
  };
}

#endif // DUNE_BASISMATRIX_HH
