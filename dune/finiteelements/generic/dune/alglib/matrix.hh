// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALGLIB_MATRIX_HH
#define DUNE_ALGLIB_MATRIX_HH

#include <cassert>

#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/vector.hh>
#include <dune/common/field.hh>

#include <alglib/ap.h>
#include <alglib/inv.h>
// #include <alglib/rcond.h>
// #include <alglib/sevd.h>

namespace Dune
{

  namespace AlgLib
  {

    template< class F, bool aligned = false >
    class Matrix;


    template< unsigned int precision, bool aligned >
    class Matrix< MultiPrecision< precision >, aligned >
    {
      typedef Matrix< MultiPrecision< precision >, aligned > This;

    public:
      typedef MultiPrecision< precision > Field;
      typedef AlgLib::Vector< Field > Vector;

    private:
      typedef amp::ampf< precision > RealField;
      typedef ap::template_2d_array< RealField, aligned > RealMatrix;

    public:
      operator const RealMatrix & () const
      {
        return matrix_;
      }

      operator RealMatrix & ()
      {
        return matrix_;
      }

      const RealMatrix &amp () const
      {
        return matrix_;
      }

      const Field &operator() ( const unsigned int row, const unsigned int col ) const
      {
        return static_cast< const Field & >( matrix_( row, col ) );
      }

      Field &operator() ( const unsigned int row, const unsigned int col )
      {
        assert(row<rows());
        assert(col<cols());
        return static_cast< Field & >( matrix_( row, col ) );
      }

      unsigned int rows () const
      {
        return matrix_.gethighbound( 1 )+1;
      }

      unsigned int cols () const
      {
        return matrix_.gethighbound( 2 )+1;
      }

      const Field *rowPtr ( const unsigned int row ) const
      {
        assert(row<rows());
        const int lastCol = matrix_.gethighbound( 2 );
        ap::const_raw_vector< RealField > rowVector = matrix_.getrow( row, 0, lastCol );
        assert( (rowVector.GetStep() == 1) && (rowVector.GetLength() == lastCol+1) );
        return static_cast< const Field * >( rowVector.GetData() );
      }

      Field *rowPtr ( const unsigned int row )
      {
        assert(row<rows());
        const int lastCol = matrix_.gethighbound( 2 );
        ap::raw_vector< RealField > rowVector = matrix_.getrow( row, 0, lastCol );
        assert( (rowVector.GetStep() == 1) && (rowVector.GetLength() == lastCol+1) );
        return static_cast< Field * >( rowVector.GetData() );
      }

      void resize ( const unsigned int rows, const unsigned int cols )
      {
        matrix_.setbounds( 0, rows-1, 0, cols-1 );
      }

      bool invert ()
      {
        assert( rows() == cols() );
        return inv::rmatrixinverse< precision >( matrix_, rows() );
      }

#if 0
      Field conditionOne () const
      {
        return Field( 1 ) / rcond::rmatrixrcond1( matrix_, rows() );
      }

      Field conditionTwo () const
      {
        const unsigned int n = rows();
        Vector d( n );
        This z;
        if( !sevd::smatrixevd< precision >( (*this), n, 0, false, d, z ) )
          DUNE_THROW( MathError, "Eigenvalue computation failed." );
        return (d[ n-1 ] / d[ 0 ]);
      }

      Field conditionInfty () const
      {
        return Field( 1 ) / rcond::rmatrixrcondinf( matrix_, rows() );
      }

      // note that the sparse matrix is assumed quadratic, here
      template< class SparseMatrix >
      void fillFromSparseMatrix ( const SparseMatrix &sparseMatrix )
      {
        const unsigned int rows = sparseMatrix.rows();
        resize( rows, rows );
        for( unsigned int i = 0; i < rows; ++i )
        {
          for( unsigned int j = 0; j < rows; ++j )
            (*this)( i, j ) = Field( 0 );

          const unsigned int nonZero = sparseMatrix.nonZero( i );
          for( unsigned int k = 0; k < nonZero; ++k )
          {
            if( !sparseMatrix.isZero( i, k ) )
            {
              const unsigned int j = sparseMatrix.column( i, k );
              field_cast( sparseMatrix( i, k ), (*this)( i, j ) );
            }
          }
        }
      }
#endif

    private:
      RealMatrix matrix_;
    };



    template< class Field, bool aligned >
    inline std::ostream &operator<<(std::ostream &out, const Matrix<Field,aligned> &mat)
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

}

#endif
