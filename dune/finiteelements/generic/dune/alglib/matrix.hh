// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALGLIB_MATRIX_HH
#define DUNE_ALGLIB_MATRIX_HH

#include <cassert>

#include <dune/alglib/multiprecision.hh>

#include <alglib/ap.h>
#include <alglib/inv.h>

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

      const Field &operator() ( const unsigned int row, const unsigned int col ) const
      {
        return static_cast< const Field & >( matrix_( row, col ) );
      }

      Field &operator() ( const unsigned int row, const unsigned int col )
      {
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
        const int lastCol = matrix_.gethighbound( 2 );
        ap::const_raw_vector< RealField > rowVector = matrix_.getrow( row, 0, lastCol );
        assert( (rowVector.GetStep() == 1) && (rowVector.GetLength() == cols()) );
        return static_cast< const Field * >( rowVector.GetData() );
      }

      Field *rowPtr ( const unsigned int row )
      {
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
        assert( matrix_.rows() == matrix_.cols() );
        return inv::rmatrixinverse< precision >( matrix_, matrix_.rows() );
      }

    private:
      RealMatrix matrix_;
    };

  }

}

#endif
