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
  template< class V >
  class CoeffMatrix
  {
    typedef V Vector;
    typedef typename Vector::field_type Field;

  public:
    CoeffMatrix()
      : coeff_(0),
        rows_(0),
        numRows_(0)
    {}

    ~CoeffMatrix()
    {
      delete [] coeff_;
      delete [] rows_;
    }

    const int size () const
    {
      return numRows_;
    }

    template <class RangeVector>
    void print(std::ostream& out,
               std::vector< RangeVector > &x) const {
      /*
         size_t numLsg = numRows_;
         Field *row = rows_[0];
         for( unsigned int r=0;r<numLsg;++r )
         {
         out << "f_" << r << "(" << char('a');
         for (int i=1;i<RangeVector::field_type::dimension;++i)
          out << "," << char('a'+i);
         out << ")=";
         RangeVector *itx = (&x[0]);
         bool first = true;
         for (; row != rows_[r+1]; ++row, ++itx) {
          if (*row > 1e-15) {
            out << ((!first)?" + ":"") << (*row) << "*" << (*itx);
            first = false;
          }
          else if (*row < -1e-15) {
            out << " - " << -(*row) << "*" << (*itx);
            first = false;
          }
         }
         out << std::endl;
         }
       */
    }

    template< class DomainVector, class RangeVector >
    void mult ( const std::vector< DomainVector > &x, std::vector< RangeVector > &y ) const
    {
      size_t numLsg = y.size();
      assert( numLsg <= (size_t)numRows_ );
      Vector *row = rows_[ 0 ];
      for( size_t r = 0; r < numLsg; ++r )
      {
        Vector val = Zero<Vector>();
        const DomainVector *itx = &(x[ 0 ]);
        for( ; row != rows_[ r+1 ]; ++row, ++itx)
          val += (*row) * (*itx) ;
        field_cast(val,y[r]);
      }
    }

    template< class FullMatrix >
    void fill ( const FullMatrix &mat )
    {
      numRows_ = mat.rowSize();
      int size = 0;
      for( int r = 0; r < numRows_; ++r )
        size += mat.colSize( r );

      coeff_ = new Vector[ size ];
      rows_ = new Vector*[ numRows_+1 ];
      rows_[ 0 ] = coeff_;
      for( int r = 0; r < numRows_; ++r )
      {
        rows_[ r+1 ] = rows_[ r ] + mat.colSize( r );
        int c = 0;
        for( Vector *it = rows_[ r ]; it != rows_[ r+1 ]; ++it, ++c )
          field_cast(mat(r,c),*it);
      }
    }

  private:
    Vector *coeff_;
    Vector **rows_;
    int numRows_;
  };
}

#endif // DUNE_COEFFMATRIX_HH
