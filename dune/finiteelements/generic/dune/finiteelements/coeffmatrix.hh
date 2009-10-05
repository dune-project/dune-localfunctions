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
  template <class Vector1,class Vector2>
  struct Mult;
  template <class Field,int dimRange,int dimDomain>
  struct Mult<FieldVector<Field,dimRange>,FieldVector<Field,dimDomain> >
  {
    typedef FieldVector<Field,dimRange> Vector1;
    typedef FieldVector<Field,dimDomain> Vector2;
    typedef FieldMatrix<Field,dimRange,dimRange> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      for (int r=0; r<dimRange; ++r)
        for (int d=0; d<dimDomain; ++d)
          res[r][d] += vec1[r]*vec2[d];
    }
  };
  template <class Field,int dimRange>
  struct Mult<FieldVector<Field,dimRange>,Field>
  {
    typedef FieldVector<Field,dimRange> Vector1;
    typedef Field Vector2;
    typedef FieldVector<Field,dimRange> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      for (int r=0; r<dimRange; ++r)
        res[r] += vec1[r]*vec2;
    }
  };


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
    void mult ( const std::vector< DomainVector > &x,
                std::vector< RangeVector > &y ) const
    {
      typedef Mult<Vector,DomainVector> Multiply;
      typedef typename Multiply::Result Result;
      size_t numLsg = y.size();
      assert( numLsg <= (size_t)numRows_ );
      Vector *row = rows_[ 0 ];
      for( size_t r = 0; r < numLsg; ++r )
      {
        Result val(0.);
        const DomainVector *itx = &(x[ 0 ]);
        int i=0;
        for( ; row != rows_[ r+1 ]; ++row, ++itx, ++i)
        {
          assert(i<x.size());
          Multiply::add(*row,*itx,val);
          // val += (*itx) * (*row) ;
        }
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
