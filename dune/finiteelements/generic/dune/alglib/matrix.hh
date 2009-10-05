// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALGLIB_MATRIX_HH
#define DUNE_ALGLIB_MATRIX_HH

#include <alglib/ap.h>
#include <alglib/amp.h>
#include <alglib/inv.h>

namespace Dune
{

  namespace AlgLib
  {

    template< unsigned int precision, bool aligned = false >
    class Matrix
    {
      typedef Matrix< precision, aligned > This;

    public:
      //typedef amp::ampf< precision > Field;
      typedef double Field;

    private:
      typedef ap::template_2d_array< Field, aligned > RealMatrix;

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
        return matrix_( row, col );
      }

      Field &operator() ( const unsigned int row, const unsigned int col )
      {
        return matrix_( row, col );
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
        ap::const_raw_vector< Field > rowVector = matrix_.getrow( row, 0, lastCol );
        assert( (rowVector.GetStep() == 1) && (rowVector.GetLength() == cols()) );
        return rowVector.GetData();
      }

      Field *rowPtr ( const unsigned int row )
      {
        const int lastCol = matrix_.gethighbound( 2 );
        ap::raw_vector< Field > rowVector = matrix_.getrow( row, 0, lastCol );
        assert( (rowVector.GetStep() == 1) && (rowVector.GetLength() == cols()) );
        return rowVector.GetData();
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
