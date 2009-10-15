// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIPRECISION_HH
#define DUNE_MULTIPRECISION_HH

#include <alglib/amp.h>

#if 0
namespace Dune
{

  namespace AlgLib
  {

    template< unsigned int precision >
    class MultiPrecision
      : public amp::ampf< precision >
    {
      typedef amp::ampf< precision > Base;

    public:
      MultiPrecision ()
        : Base()
      {}

      template< class T >
      MultiPrecision ( const T &v )
        : Base( v )
      {}
    };



    template< unsigned int precision >
    inline MultiPrecision< precision >
    operator+ ( const MultiPrecision< precision > &a, const MultiPrecision< precision > &b )
    {
      typedef amp::ampf< precision > F;
      return ((const F &)a + (const F &)b);
    }

    template< unsigned int precision >
    inline MultiPrecision< precision >
    operator- ( const MultiPrecision< precision > &a, const MultiPrecision< precision > &b )
    {
      typedef amp::ampf< precision > F;
      return ((const F &)a - (const F &)b);
    }

    template< unsigned int precision >
    inline MultiPrecision< precision >
    operator- ( const MultiPrecision< precision > &a )
    {
      typedef amp::ampf< precision > F;
      return -((const F &)a);
    }

    template< unsigned int precision >
    inline MultiPrecision< precision >
    operator* ( const MultiPrecision< precision > &a, const MultiPrecision< precision > &b )
    {
      typedef amp::ampf< precision > F;
      return ((const F &)a * (const F &)b);
    }

    template< unsigned int precision >
    inline MultiPrecision< precision >
    operator/ ( const MultiPrecision< precision > &a, const MultiPrecision< precision > &b )
    {
      typedef amp::ampf< precision > F;
      return ((const F &)a / (const F &)b);
    }



    template< unsigned int precision >
    inline std::ostream &
    operator<< ( std::ostream &out, const MultiPrecision< precision > &value )
    {
      return out << value.toDec();
    }

  }

}



namespace std
{

  template< unsigned int precision >
  inline Dune::AlgLib::MultiPrecision< precision >
  sqrt ( const Dune::AlgLib::MultiPrecision< precision > &a )
  {
    return amp::sqrt( a );
  }

}
#endif

namespace std
{

  template< unsigned int precision >
  inline std::ostream &
  operator<< ( std::ostream &out, const amp::ampf< precision > &value )
  {
    return out << value.toDec();
  }

  template< unsigned int precision >
  inline amp::ampf< precision > sqrt ( const amp::ampf< precision > &a )
  {
    return amp::sqrt( a );
  }

}


#endif // #ifndef DUNE_MULTIPRECISION_HH
