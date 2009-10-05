// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALGLIB_VECTOR_HH
#define DUNE_ALGLIB_VECTOR_HH

#include <cassert>

#include <dune/alglib/multiprecision.hh>

#include <alglib/ap.h>

namespace Dune
{

  namespace AlgLib
  {

    template< class F, bool aligned = false >
    class Vector;


    template< unsigned int precision, bool aligned >
    class Vector< MultiPrecision< precision >, aligned >
    {
      typedef Vector< MultiPrecision< precision >, aligned > This;

    public:
      typedef MultiPrecision< precision > Field;

    private:
      typedef amp::ampf< precision > RealField;
      typedef ap::template_1d_array< RealField, aligned > RealVector;

    public:
      Vector ()
      {}

      Vector ( unsigned int size )
      {
        resize( size );
      }

      operator const RealVector & () const
      {
        return vector_;
      }

      operator RealVector & ()
      {
        return vector_;
      }

      const Field &operator[] ( const unsigned int i ) const
      {
        return static_cast< const Field & >( vector_( i ) );
      }

      Field &operator[] ( const unsigned int i )
      {
        return static_cast< Field & >( vector_( i ) );
      }

      unsigned int size () const
      {
        return vector_.gethighbound()+1;
      }

      const Field *ptr ( const unsigned int row ) const
      {
        const int lastIndex = vector_.gethighbound();
        ap::const_raw_vector< RealField > rawVector = vector_.getvector( 0, lastIndex );
        assert( (rawVector.GetStep() == 1) && (rawVector.GetLength() == lastIndex+1) );
        return static_cast< const Field * >( rawVector.GetData() );
      }

      Field *ptr ( const unsigned int row )
      {
        const int lastIndex = vector_.gethighbound();
        ap::raw_vector< RealField > rawVector = vector_.getvector( 0, lastIndex );
        assert( (rawVector.GetStep() == 1) && (rawVector.GetLength() == lastIndex+1) );
        return static_cast< Field * >( rawVector.GetData() );
      }

      void resize ( const unsigned int size )
      {
        vector_.setbounds( 0, size-1 );
      }

    private:
      RealVector vector_;
    };

  }

}

#endif // #ifndef DUNE_ALGLIB_VECTOR_HH
