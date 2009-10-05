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
  template <class MatrixEntry,class BasisEntry>
  struct Mult;

  template <class Field,int dimRange>
  struct Mult< Field,FieldVector<Field,dimRange> >
  {
    typedef FieldVector<Field,dimRange> BasisEntry;
    static void add(const Field &vec1, const BasisEntry &vec2,
                    BasisEntry &res)
    {
      res.axpy(vec1,vec2);
    }
  };
  template <class Field>
  struct Mult< Field,Field >
  {
    typedef Field BasisEntry;
    static void add(const Field &vec1, const BasisEntry &vec2,
                    BasisEntry &res)
    {
      res += vec1*vec2;
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

    template <class RangeVector>
    void print(std::ostream& out,
               std::vector< RangeVector > &x,
               unsigned int numLsg = size()) const {
      /*
         Vector *row = rows_[0];
         for( unsigned int r=0;r<numLsg;++r )
         {
          out << "f_" << r << "(" << char('a');
          for (int i=1;i<RangeVector::dimension;++i)
            out << "," << char('a'+i);
          out << ")=";
          RangeVector *itx = (&x[0]);
          bool first = true;
          for (; row != rows_[r+1]; ++row, ++itx) {
            if ((*row)[0][0] > 1e-15) {
              out << ((!first)?" + ":"") << (*row)[0][0] << "*" << (*itx);
              first = false;
            }
            else if ((*row)[0][0] < -1e-15) {
              out << " - " << -((*row)[0][0]) << "*" << (*itx);
              first = false;
            }
          }
          out << std::endl;
         }
       */
    }

    template< class BasisVector, class RangeVector >
    void mult ( const BasisVector &x,
                std::vector< RangeVector > &y ) const
    {
      typedef typename BasisVector::RangeVector DomainVector;
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
          Multiply::add(*row,*itx,val);
        }
        field_cast(val,y[r]);
      }
    }

    template< class FullMatrix >
    void fill ( const FullMatrix &mat )
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

      coeff_ = new Field[size ];
      rows_ = new Field*[ numRows_+1 ];
      rows_[ 0 ] = coeff_;
      for( int r = 0; r < numRows_; ++r )
      {
        rows_[ r+1 ] = rows_[ r ] + mat.colSize( r );
        int c = 0;
        for( Field *it = rows_[ r ]; it != rows_[ r+1 ]; ++it, ++c )
          field_cast(mat(r,c),*it);
      }
    }

  private:
    CoeffMatrix(const CoeffMatrix&);
    CoeffMatrix &operator= (const CoeffMatrix&);
    Field *coeff_;
    Field **rows_;
    int numRows_,numCols_;
  };

}

#endif // DUNE_COEFFMATRIX_HH
