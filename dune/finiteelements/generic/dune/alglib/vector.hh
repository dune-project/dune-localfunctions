// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALGLIB_VECTOR_HH
#define DUNE_ALGLIB_VECTOR_HH

#include <cassert>
#include <vector>

#if HAVE_ALGLIB
#include <alglib/amp.h>
#endif

namespace Dune
{

  namespace AlgLib
  {

    template< class F, bool aligned = false >
    class Vector;

    template< class F, bool aligned >
    class Vector
    {
      typedef Vector< F, aligned > This;

    public:
      typedef F Field;

    private:
      typedef std::vector< Field > RealVector;

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
        return vector_[ i ];
      }

      Field &operator[] ( const unsigned int i )
      {
        return vector_[ i ];
      }

      unsigned int size () const
      {
        return vector_.size();
      }

      const Field *ptr ( const unsigned int row ) const
      {
        return &(vector_[row]);
      }

      Field *ptr ( const unsigned int row )
      {
        return &(vector_[row]);
      }

      void resize ( const unsigned int size )
      {
        vector_.resize(size);
      }

    private:
      RealVector vector_;
    };

#if HAVE_ALGLIB
    template< unsigned int precision, bool aligned >
    class Vector< amp::ampf< precision >, aligned >
    {
      typedef Vector< amp::ampf< precision >, aligned > This;

    public:
      typedef amp::ampf< precision > Field;

    private:
      typedef ap::template_1d_array< Field, aligned > RealVector;

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
        return vector_( i );
      }

      Field &operator[] ( const unsigned int i )
      {
        return vector_( i );
      }

      unsigned int size () const
      {
        return vector_.gethighbound()+1;
      }

      const Field *ptr ( const unsigned int row ) const
      {
        const int lastIndex = vector_.gethighbound();
        ap::const_raw_vector< Field > rawVector = vector_.getvector( 0, lastIndex );
        assert( (rawVector.GetStep() == 1) && (rawVector.GetLength() == lastIndex+1) );
        return rawVector.GetData();
      }

      Field *ptr ( const unsigned int row )
      {
        const int lastIndex = vector_.gethighbound();
        ap::raw_vector< Field > rawVector = vector_.getvector( 0, lastIndex );
        assert( (rawVector.GetStep() == 1) && (rawVector.GetLength() == lastIndex+1) );
        return rawVector.GetData();
      }

      void resize ( const unsigned int size )
      {
        vector_.setbounds( 0, size-1 );
      }

    private:
      RealVector vector_;
    };
#endif
  }

}

#endif // #ifndef DUNE_ALGLIB_VECTOR_HH
