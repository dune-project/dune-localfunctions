// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIPRECISION_HH
#define DUNE_MULTIPRECISION_HH

#include <alglib/amp.h>

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

  }

}

#endif // #ifndef DUNE_MULTIPRECISION_HH
