// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_COEFFMATRIX_HH
#define DUNE_COEFFMATRIX_HH
#include <cassert>
#include <iostream>
#include <vector>
#include <dune/common/fvector.hh>
#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/tensor.hh>

namespace Dune
{
  /*************************************************
  * Default class for storing a coefficient matrix
  * for the PolynomialBasis. Basically a simple
  * CRS structure is used. The additional complexity
  * is due to the storage and efficient evaluation
  * of higher order derivatives. See the remarks
  * in tensor.hh which also hold true for this file.
  *************************************************/
  template <class Field, class Field2>
  struct Mult
  {
    typedef Field2 BasisEntry;
    static void add(const Field &vec1, const BasisEntry &vec2,
                    BasisEntry &res)
    {
      res += vec1*vec2;
    }
  };

  template <class Field,class Field2, int dimRange>
  struct Mult< Field,FieldVector<Field2,dimRange> >
  {
    typedef FieldVector<Field2,dimRange> BasisEntry;
    static void add(const Field &vec1, const BasisEntry &vec2,
                    BasisEntry &res)
    {
      res.axpy(vec1,vec2);
    }
  };

  template< class F , unsigned int bSize >
  class SparseCoeffMatrix
  {
  public:
    typedef F Field;
    static const unsigned int blockSize = bSize;
    typedef SparseCoeffMatrix<Field,blockSize> This;

    SparseCoeffMatrix()
      : coeff_(0),
        rows_(0),
        skip_(0),
        numRows_(0),
        numCols_(0)
    {}

    ~SparseCoeffMatrix()
    {
      delete [] coeff_;
      delete [] rows_;
      delete [] skip_;
    }

    unsigned int size () const
    {
      return numRows_/blockSize;
    }
    unsigned int baseSize () const
    {
      return numCols_;
    }

    template< class BasisIterator, class FF>
    void mult ( const BasisIterator &x,
                unsigned int numLsg,
                FF *y ) const
    {
      typedef typename BasisIterator::Derivatives XDerivatives;
      assert( numLsg*blockSize <= (size_t)numRows_ );
      unsigned int row = 0;
      Field *pos = rows_[ 0 ];
      unsigned int *skipIt = skip_;
      XDerivatives val;
      for( size_t i = 0; i < numLsg; ++i)
      {
        for( unsigned int r = 0; r < blockSize; ++r, ++row )
        {
          val = 0;
          BasisIterator itx = x;
          for( ; pos != rows_[ row+1 ]; ++pos, ++skipIt )
          {
            itx += *skipIt;
            val.axpy(*pos,*itx);
          }
          DerivativeAssign<XDerivatives,FF>::apply(r,val,*(y+i*XDerivatives::size*blockSize));
        }
      }
    }
    template< class BasisIterator, class Vector>
    void mult ( const BasisIterator &x,
                Vector &y ) const
    {
      typedef typename Vector::value_type YDerivatives;
      typedef typename BasisIterator::Derivatives XDerivatives;
      size_t numLsg = y.size();
      assert( numLsg*blockSize <= (size_t)numRows_ );
      unsigned int row = 0;
      Field *pos = rows_[ 0 ];
      unsigned int *skipIt = skip_;
      XDerivatives val;
      for( size_t i = 0; i < numLsg; ++i)
      {
        for( unsigned int r = 0; r < blockSize; ++r, ++row )
        {
          val = 0;
          BasisIterator itx = x;
          for( ; pos != rows_[ row+1 ]; ++pos, ++skipIt )
          {
            itx += *skipIt;
            val.axpy(*pos,*itx);
          }
          DerivativeAssign<XDerivatives,YDerivatives>::apply(r,val,y[i]);
        }
      }
    }
    template <unsigned int deriv, class BasisIterator, class Vector>
    void mult ( const BasisIterator &x,
                Vector &y ) const
    {
      typedef typename Vector::value_type YDerivatives;
      typedef typename BasisIterator::Derivatives XDerivatives;
      typedef FieldVector<typename XDerivatives::Field,YDerivatives::dimension> XLFETensor;
      size_t numLsg = y.size();
      assert( numLsg*blockSize <= (size_t)numRows_ );
      unsigned int row = 0;
      Field *pos = rows_[ 0 ];
      unsigned int *skipIt = skip_;
      for( size_t i = 0; i < numLsg; ++i)
      {
        XLFETensor val(typename XDerivatives::Field(0));
        for( unsigned int r = 0; r < blockSize; ++r, ++row )
        {
          BasisIterator itx = x;
          for( ; pos != rows_[ row+1 ]; ++pos, ++skipIt )
          {
            itx += *skipIt;
            LFETensorAxpy<XDerivatives,XLFETensor,deriv>::apply(r,*pos,*itx,val);
          }
        }
        field_cast(val,y[i]);
      }
    }

    template< class RowMatrix >
    void fill ( const RowMatrix &mat, bool verbose=false )
    {
      numRows_ = mat.rows();
      numCols_ = mat.cols();
      unsigned int size = numRows_*numCols_;

      delete [] coeff_;
      delete [] rows_;
      delete [] skip_;

      Field* coeff = new Field[ size ];
      // we always initialize the next skip entry to zero,
      // including the one following the end, so allocate
      // size+1 entries so we will stay within the bounds.
      unsigned int *skip = new unsigned int[ size+1 ];
      rows_ = new Field*[ numRows_+1 ];
      std::vector<Field> row( numCols_ );

      rows_[ 0 ] = coeff;
      Field *cit = coeff;
      unsigned int *sit = skip;
      for( unsigned int r = 0; r < numRows_; ++r )
      {
        *sit = 0;
        mat.row( r, row );
        for( unsigned int c = 0; c < numCols_; ++c )
        {
          const Field &val = row[c];
          if (val < Zero<Field>() || Zero<Field>() < val)
          {
            *cit = val;
            ++sit;
            ++cit;
            *sit = 1;
          } else
          {
            ++(*sit);
          }
        }
        rows_[ r+1 ] = cit;
      }
      assert( size_t(rows_[numRows_]-rows_[0]) <= size_t(size) );
      size = rows_[numRows_]-rows_[0];
      coeff_ = new Field[ size ];
      skip_ = new unsigned int[ size ];
      for (unsigned int i=0; i<size; ++i)
      {
        coeff_[i] = coeff[i];
        skip_[i] = skip[i];
      }
      for (unsigned int i=0; i<=numRows_; ++i)
        rows_[ i ] = coeff_ + (rows_[ i ] - coeff);

      delete [] coeff;
      delete [] skip;

      if (verbose)
        std::cout << "Entries: " << (rows_[numRows_]-rows_[0])
                  << " full: " << numCols_*numRows_
                  << std::endl;
    }
    // b += a*C[k]
    template <class Vector>
    void addRow( unsigned int k, const Field &a, Vector &b) const
    {
      assert(k<numRows_);
      unsigned int j=0;
      unsigned int *skipIt = skip_ + (rows_[ k ]-rows_[ 0 ]);
      for( Field *pos = rows_[ k ];
           pos != rows_[ k+1 ];
           ++pos, ++skipIt )
      {
        j += *skipIt;
        assert( j < b.size() );
        b[j] += field_cast<typename Vector::value_type>( (*pos)*a );  // field_cast
      }
    }
  private:
    SparseCoeffMatrix ( const This &other )
      : numRows_( other.numRows_ ),
        numCols_( other.numCols_ )
    {
      const unsigned int size = other.rows_[numRows_]-other.rows_[0];
      coeff_ = new Field[ size ];
      rows_ = new Field*[ numRows_+1 ];
      skip_ = new unsigned int[ size ];
      for (unsigned int i=0; i<size; ++i)
      {
        coeff_[i] = other.coeff_[i];
        skip_[i] = other.skip_[i];
      }
      for (unsigned int i=0; i<=numRows_; ++i)
        rows_[ i ] = coeff_ + (other.rows_[ i ] - other.coeff_);
    }

    This &operator= (const This&);
    Field *coeff_;
    Field **rows_;
    unsigned int *skip_;
    unsigned int numRows_,numCols_;
  };

}

#endif // DUNE_COEFFMATRIX_HH
