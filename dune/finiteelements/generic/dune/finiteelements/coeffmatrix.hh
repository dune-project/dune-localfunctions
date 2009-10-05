// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_COEFFMATRIX_HH
#define DUNE_COEFFMATRIX_HH
#include <cassert>
#include <iostream>
#include <vector>
#include <dune/common/fvector.hh>
#include <dune/common/field.hh>

namespace Dune
{
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

  template< class F >
  class CoeffMatrix
  {
  public:
    typedef CoeffMatrix<F> This;
    typedef F Field;

    CoeffMatrix()
      : coeff_(0),
        rows_(0),
        numRows_(0),
        numCols_(0)
    {}

    ~CoeffMatrix()
    {
      delete [] coeff_;
      delete [] rows_;
    }

    const unsigned int size () const
    {
      return numRows_;
    }
    const unsigned int baseSize () const
    {
      return numCols_;
    }

    template< class BasisVector, class RangeVector >
    void mult ( const BasisVector &x,
                std::vector< RangeVector > &y ) const
    {
      typedef typename BasisVector::Block DomainVector;
      typedef Mult<Field,DomainVector> Multiply;
      size_t numLsg = y.size();
      assert( numLsg <= (size_t)numRows_ );
      Field *row = rows_[ 0 ];
      for( size_t r = 0; r < numLsg; ++r )
      {
        DomainVector val(0.);
        // typename BasisVector::Iterator itx = x.begin();
        BasisVector itx = x;
        for( ; row != rows_[ r+1 ]; ++row, ++itx )
        {
          Multiply::add(*row,itx.block(),val);
        }
        field_cast(val,y[r]);
      }
    }

    template< class FullMatrix >
    void fill ( const FullMatrix &mat, bool verbose = false )
    {
      unsigned int zeros = 0;
      numRows_ = mat.rowSize();
      numCols_ = 0;
      int size = 0;
      for( int r = 0; r < numRows_; ++r ) {
        size += mat.colSize( r );
        numCols_ = std::max( numCols_ , mat.colSize( r ) );
      }

      delete [] coeff_;
      delete [] rows_;

      coeff_ = new Field[size ];
      rows_ = new Field*[ numRows_+1 ];
      rows_[ 0 ] = coeff_;
      for( int r = 0; r < numRows_; ++r )
      {
        rows_[ r+1 ] = rows_[ r ] + mat.colSize( r );
        int c = 0;
        for( Field *it = rows_[ r ]; it != rows_[ r+1 ]; ++it, ++c )
        {
          Field val;
          field_cast(mat(r,c),val);
          if (std::abs(val)<1e-10)
            ++zeros;
          *it = val;
        }
      }
      if (verbose)
        std::cout << "Entries: " << (rows_[numRows_]-rows_[0]) << " "
                  << "Zeros: " << zeros << " "
                  << "Precentage: " << 100.*double(zeros)/double(rows_[numRows_]-rows_[0])
                  << std::endl;
    }

  private:
    CoeffMatrix(const CoeffMatrix&);
    CoeffMatrix &operator= (const CoeffMatrix&);
    Field *coeff_;
    Field **rows_;
    int numRows_,numCols_;
  };

  template< class F >
  class SparseCoeffMatrix
  {
  public:
    typedef SparseCoeffMatrix<F> This;
    typedef F Field;

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

    const unsigned int size () const
    {
      return numRows_;
    }
    const unsigned int baseSize () const
    {
      return numCols_;
    }

    template< class BasisVector, class Vector >
    void mult ( const BasisVector &x,
                Vector  &y ) const
    {
      typedef typename BasisVector::Block DomainVector;
      typedef Mult<Field,DomainVector> Multiply;
      size_t numLsg = y.size();
      assert( numLsg <= (size_t)numRows_ );
      Field *row = rows_[ 0 ];
      unsigned int *skipIt = skip_;
      for( size_t r = 0; r < numLsg; ++r )
      {
        DomainVector val(0.);
        // typename BasisVector::Iterator itx = x.begin();
        BasisVector itx = x;
        for( ; row != rows_[ r+1 ]; ++row, ++skipIt )
        {
          itx += *skipIt;
          Multiply::add(*row,itx.block(),val);
        }
        field_cast(val,y[r]);
      }
    }

    template< class FullMatrix >
    void fill ( const FullMatrix &mat, bool verbose=false )
    {
      numRows_ = mat.rowSize();
      numCols_ = 0;
      int size = 0;
      for( int r = 0; r < numRows_; ++r ) {
        size += mat.colSize( r );
        numCols_ = std::max( numCols_ , mat.colSize( r ) );
      }

      delete [] coeff_;
      delete [] rows_;
      delete [] skip_;

      coeff_ = new Field[ size ];
      rows_ = new Field*[ numRows_+1 ];
      skip_ = new unsigned int[ size+1 ];
      rows_[ 0 ] = coeff_;
      Field *cit = coeff_;
      unsigned int *sit = skip_;
      for( int r = 0; r < numRows_; ++r )
      {
        *sit = 0;
        for( int c = 0; c < mat.colSize( r ); ++c )
        {
          Field val;
          field_cast(mat(r,c),val);
          if (std::abs(val)<1e-10)
          {
            ++(*sit);
          } else
          {
            *cit = val;
            ++sit;
            ++cit;
            *sit = 1;
          }
        }
        rows_[ r+1 ] = cit;
      }
      if (verbose)
        std::cout << "Entries: " << (rows_[numRows_]-rows_[0]) << std::endl;
    }

  private:
    SparseCoeffMatrix(const This&);
    This &operator= (const This&);
    Field *coeff_;
    Field **rows_;
    unsigned int *skip_;
    int numRows_,numCols_;
  };

}

#endif // DUNE_COEFFMATRIX_HH
