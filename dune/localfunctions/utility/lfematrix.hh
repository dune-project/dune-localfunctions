// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_UTILITY_LFEMATRIX_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_LFEMATRIX_HH

#include <cassert>
#include <vector>

#include "field.hh"

namespace Dune
{

  template< class F >
  class LFEMatrix
  {
    typedef LFEMatrix< F > This;
    typedef std::vector< F > Row;
    typedef std::vector<Row> RealMatrix;

  public:
    typedef F Field;

    operator const RealMatrix & () const
    {
      return matrix_;
    }

    operator RealMatrix & ()
    {
      return matrix_;
    }

    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      assert(row<rows());
      for (int i=0; i<cols(); ++i)
        field_cast(matrix_[row][i], vec[i]);
    }

    const Field &operator() ( const unsigned int row, const unsigned int col ) const
    {
      assert(row<rows());
      assert(col<cols());
      return matrix_[ row ][ col ];
    }

    Field &operator() ( const unsigned int row, const unsigned int col )
    {
      assert(row<rows());
      assert(col<cols());
      return matrix_[ row ][ col ];
    }

    unsigned int rows () const
    {
      return rows_;
    }

    unsigned int cols () const
    {
      return cols_;
    }

    const Field *rowPtr ( const unsigned int row ) const
    {
      assert(row<rows());
      return &(matrix_[row][0]);
    }

    Field *rowPtr ( const unsigned int row )
    {
      assert(row<rows());
      return &(matrix_[row][0]);
    }

    void resize ( const unsigned int rows, const unsigned int cols )
    {
      matrix_.resize(rows);
      for (unsigned int i=0; i<rows; ++i)
        matrix_[i].resize(cols);
      rows_ = rows;
      cols_ = cols;
    }

    bool invert ()
    {
      using std::abs;
      assert( rows() == cols() );
      std::vector<unsigned int> p(rows());
      for (unsigned int j=0; j<rows(); ++j)
        p[j] = j;
      for (unsigned int j=0; j<rows(); ++j)
      {
        // pivot search
        unsigned int r = j;
        Field max = abs( (*this)(j,j) );
        for (unsigned int i=j+1; i<rows(); ++i)
        {
          if ( abs( (*this)(i,j) ) > max )
          {
            max = abs( (*this)(i,j) );
            r = i;
          }
        }
        if (max == Zero<Field>())
          return false;
        // row swap
        if (r > j)
        {
          for (unsigned int k=0; k<cols(); ++k)
            std::swap( (*this)(j,k), (*this)(r,k) );
          std::swap( p[j], p[r] );
        }
        // transformation
        Field hr = Unity<Field>()/(*this)(j,j);
        for (unsigned int i=0; i<rows(); ++i)
          (*this)(i,j) *= hr;
        (*this)(j,j) = hr;
        for (unsigned int k=0; k<cols(); ++k)
        {
          if (k==j) continue;
          for (unsigned int i=0; i<rows(); ++i)
          {
            if (i==j) continue;
            (*this)(i,k) -= (*this)(i,j)*(*this)(j,k);
          }
          (*this)(j,k) *= -hr;
        }
      }
      // column exchange
      Row hv(rows());
      for (unsigned int i=0; i<rows(); ++i)
      {
        for (unsigned int k=0; k<rows(); ++k)
          hv[ p[k] ] = (*this)(i,k);
        for (unsigned int k=0; k<rows(); ++k)
          (*this)(i,k) = hv[k];
      }
      return true;
    }

  private:
    RealMatrix matrix_;
    unsigned int cols_,rows_;
  };

  template< class Field >
  inline std::ostream &operator<<(std::ostream &out, const LFEMatrix<Field> &mat)
  {
    for (unsigned int r=0; r<mat.rows(); ++r)
    {
      out << field_cast<double>(mat(r,0));
      for (unsigned int c=1; c<mat.cols(); ++c)
      {
        out << " , " << field_cast<double>(mat(r,c));
      }
      out << std::endl;
    }
    return out;
  }
}

#endif // #ifndef DUNE_LOCALFUNCTIONS_UTILITY_LFEMATRIX_HH
