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
  template <class Field,int dimRange,int dimDomain>
  struct Mult<FieldVector<Field,dimRange>,FieldVector<Field,dimDomain> >
  {
    typedef FieldVector<Field,dimRange> Vector1;
    typedef FieldVector<Field,dimDomain> Vector2;
    typedef FieldMatrix<Field,dimRange,dimDomain> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      for (int r=0; r<dimRange; ++r)
        for (int d=0; d<dimDomain; ++d)
          res[r][d] += vec1[r]*vec2[d];
    }
  };

  template <class Field>
  struct Mult<FieldMatrix<Field,1,1>,Field >
  {
    typedef FieldMatrix<Field,1,1> Vector1;
    typedef Field Vector2;
    typedef FieldVector<Field,1> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      res[0] += vec1[0][0]*vec2;
    }
  };
  template <class Field,int dimRange>
  struct Mult<FieldMatrix<Field,dimRange,dimRange>,FieldVector<Field,dimRange> >
  {
    typedef FieldMatrix<Field,dimRange,dimRange> Vector1;
    typedef FieldVector<Field,dimRange> Vector2;
    typedef FieldVector<Field,dimRange> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      for (int r=0; r<dimRange; ++r)
        for (int i=0; i<dimRange; ++i)
          res[r] += vec1[r][i]*vec2[i];
    }
  };
  template <class Field,int dimRange,int dimDomain>
  struct Mult<FieldMatrix<Field,dimRange,dimRange>,FieldMatrix<Field,dimRange,dimDomain> >
  {
    typedef FieldMatrix<Field,dimRange,dimRange> Vector1;
    typedef FieldMatrix<Field,dimRange,dimDomain> Vector2;
    typedef FieldMatrix<Field,dimRange,dimDomain> Result;
    static void add(const Vector1 &vec1, const Vector2 &vec2,
                    Result &res)
    {
      for (int r=0; r<dimRange; ++r)
        for (int d=0; d<dimDomain; ++d)
          for (int i=0; i<dimRange; ++i)
            res[r][d] += vec1[r][i]*vec2[i][d];
    }
  };


  template< class V >
  class CoeffMatrix
  {
  public:
    typedef CoeffMatrix<V> This;
    typedef V Vector;
    static const int dimension = V::dimension;
    typedef typename Vector::field_type Field;

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
      Vector *row = rows_[0];
      for( unsigned int r=0; r<numLsg; ++r )
      {
        out << "f_" << r << "(" << char('a');
        for (int i=1; i<RangeVector::dimension; ++i)
          out << "," << char('a'+i);
        out << ")=";
        RangeVector *itx = (&x[0]);
        bool first = true;
        for (; row != rows_[r+1]; ++row, ++itx) {
          if ((*row)[0][0] > 1e-15) {
            out << ((!first) ? " + " : "") << (*row)[0][0] << "*" << (*itx);
            first = false;
          }
          else if ((*row)[0][0] < -1e-15) {
            out << " - " << -((*row)[0][0]) << "*" << (*itx);
            first = false;
          }
        }
        out << std::endl;
      }
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
        for( ; row != rows_[ r+1 ]; ++row, ++itx )
        {
          Multiply::add(*row,*itx,val);
        }
        field_cast(val,y[r]);
      }
    }

    template< class BasisVector, class RangeVector >
    void mult ( const BasisVector &x,
                std::vector< RangeVector > &y ) const
    {
      typedef typename BasisVector::value_type DomainVector;
      typedef Mult<Vector,DomainVector> Multiply;
      typedef typename Multiply::Result Result;
      size_t numLsg = y.size();
      assert( numLsg <= (size_t)numRows_ );
      Vector *row = rows_[ 0 ];
      for( size_t r = 0; r < numLsg; ++r )
      {
        Result val(0.);
        typename BasisVector::const_iterator itx = x.begin();
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
    CoeffMatrix(const CoeffMatrix&);
    CoeffMatrix &operator= (const CoeffMatrix&);
    Vector *coeff_;
    Vector **rows_;
    int numRows_,numCols_;
  };

}

#endif // DUNE_COEFFMATRIX_HH
