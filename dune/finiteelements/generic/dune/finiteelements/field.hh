// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FIELD_HH
#define DUNE_FIELD_HH

namespace Dune
{

  // Unity
  // -----

  template< class Field >
  struct Unity
  {
    operator Field () const
    {
      return Field( 1 );
    }
  };



  template< class Field >
  Field operator+ ( const Unity< Field > &u, const Field &f )
  {
    return (Field)u + f;
  }

  template< class Field >
  Field operator- ( const Unity< Field > &u, const Field &f )
  {
    return (Field)u - f;
  }

  template< class Field >
  Field operator* ( const Unity< Field > &u, const Field &f )
  {
    return f;
  }

  template< class Field >
  Field operator/ ( const Unity< Field > &u, const Field &f )
  {
    return (Field)u / f;
  }



  // Zero
  // ----

  template< class Field >
  struct Zero
  {
    operator Field () const
    {
      return Field( 0 );
    }
  };

  template< class Field >
  inline bool operator< ( const Zero< Field > &, const Field &f )
  {
    return f > 1e-12;
  }

  template< class Field >
  inline bool operator< ( const Field &f, const Zero< Field > & )
  {
    return f < -1e-12;
  }

  template< class Field >
  inline bool operator> ( const Zero< Field > &z, const Field &f )
  {
    return f < z;
  }

  template< class Field >
  inline bool operator> ( const Field &f, const Zero< Field > &z )
  {
    return z < f;
  }

}

#endif // #ifndef DUNE_FIELD_HH
