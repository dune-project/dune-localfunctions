// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_UTILITY_FIELD_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_FIELD_HH

#include <dune/common/gmpfield.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  // Unity
  // -----

  /**
   * @brief A class representing the unit of a given Field
   *
   * This class can be used to assign the unit element to an
   * instance of a given Field also operators for +/- with
   * unit element are provided. Also 1/f can be evaluated.
   * Through specialization this class can be used also in the case that the
   * integer 1 is not automatically converted to the unit
   * element of the Field - the default implementation
   **/
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

  /**
   * @brief A class representing the zero of a given Field
   *
   * This class can be used to assign the zero element to an
   * instance of a given Field. An epsilon is also
   * provided for the comparison operators.
   * This class can be used also in the case that the
   * integer 0 is not automatically converted to the zero
   * element of the Field and the epsilon can be changed
   * depending on the accuracy of the Field type.
   **/
  template< class Field >
  struct Zero
  {
    operator Field () const
    {
      return Field( 0 );
    }
    static const Field epsilon()
    {
      return Field(1e-12);
    }
  };

#if HAVE_GMP
  template< unsigned int precision >
  struct Zero< GMPField< precision > >
  {
    typedef GMPField< precision > Field;
    operator Field () const
    {
      return Field( 0 );
    }
    static const Field epsilon()
    {
      return Field(1e-20);
    }
  };
#endif

  template< class Field >
  inline bool operator == ( const Zero< Field > &, const Field &f )
  {
    return ( f < Zero<Field>::epsilon() && f > -Zero<Field>::epsilon() );
  }

  template< class Field >
  inline bool operator == ( const Field &f, const Zero< Field > &z)
  {
    return ( z == f );
  }

  template< class Field >
  inline bool operator< ( const Zero< Field > &, const Field &f )
  {
    return f > Zero<Field>::epsilon();
  }

  template< class Field >
  inline bool operator< ( const Field &f, const Zero< Field > & )
  {
    return f < -Zero<Field>::epsilon();
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


  // field_cast
  // ----------

  /**
   * @brief a helper class to cast from one field
   *        to another
   *
   * This cast can be used for assignement between
   * different field types, including for example
   * between FieldVectors with different fields.
   * Specially the conversion from a special type
   * e.g. gmp to build in types are provided, the
   * other direction can be more easily handled by
   * the special field type implementation.
   **/
  template< class F2, class F1 >
  inline void field_cast ( const F1 &f1, F2 &f2 )
  {
    f2 = f1;
  }

#if HAVE_GMP
  template< unsigned int precision >
  inline void field_cast ( const Dune::GMPField< precision > &f1, double &f2 )
  {
    f2 = f1.get_d();
  }

  template< unsigned int precision >
  inline void field_cast ( const Dune::GMPField< precision > &f1, long double &f2 )
  {
    f2 = f1.get_d();
  }
#endif

  template< class F2, class F1, int dim >
  inline void field_cast ( const Dune::FieldVector< F1, dim > &f1, Dune::FieldVector< F2, dim > &f2 )
  {
    for( int d = 0; d < dim; ++d )
      field_cast( f1[ d ], f2[ d ] );
  }
  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldVector< F1, 1 > &f1, F2 &f2 )
  {
    field_cast( f1[ 0 ], f2 );
  }
  template< class F2, class F1 >
  inline void field_cast ( const F1 &f1, Dune::FieldVector< F2, 1 > &f2 )
  {
    field_cast( f1, f2[ 0 ] );
  }

  template< class F2, class F1, int rdim, int cdim >
  inline void field_cast ( const Dune::FieldMatrix< F1, rdim, cdim > &f1, Dune::FieldMatrix< F2, rdim, cdim > &f2 )
  {
    for( int r = 0; r < rdim; ++r )
      field_cast( f1[ r ], f2[ r ] );
  }
  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldMatrix<F1,1,1> &f1, Dune::FieldMatrix< F2, 1,1 > &f2 )
  {
    field_cast( f1[ 0 ][ 0 ], f2[ 0 ][ 0 ] );
  }
  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldMatrix< F1, 1,1 > &f1, F2 &f2 )
  {
    field_cast( f1[ 0 ][ 0 ], f2 );
  }
  template< class F2, class F1 >
  inline void field_cast ( const F1 &f1, Dune::FieldMatrix< F2, 1,1 > &f2 )
  {
    field_cast( f1, f2[ 0 ][ 0 ] );
  }
  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldVector<F1,1> &f1, Dune::FieldMatrix< F2, 1,1 > &f2 )
  {
    field_cast( f1[ 0 ], f2[ 0 ][ 0 ] );
  }
  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldMatrix<F1,1,1> &f1, Dune::FieldVector< F2, 1 > &f2 )
  {
    field_cast( f1[ 0 ][ 0 ], f2[ 0 ] );
  }

  template< class F2, class F1 >
  inline void field_cast ( const Dune::FieldVector< F1, 1 > &f1, Dune::FieldVector<F2, 1> &f2 )
  {
    field_cast( f1[ 0 ], f2[ 0 ] );
  }

  template< class F2,class V >
  struct FieldCast
  {
    typedef F2 type;
  };
  template< class F2,class F1,int dim >
  struct FieldCast< F2, Dune::FieldVector<F1,dim> >
  {
    typedef Dune::FieldVector<F2,dim> type;
  };
  template< class F2,class F1,int dim1, int dim2>
  struct FieldCast< F2, Dune::FieldMatrix<F1,dim1,dim2> >
  {
    typedef Dune::FieldMatrix<F2,dim1,dim2> type;
  };
  template< class F2,class V >
  inline typename FieldCast<F2,V>::type field_cast ( const V &f1 )
  {
    typename FieldCast<F2,V>::type f2;
    field_cast( f1, f2 );
    return f2;
  }


  // Precision
  // this is not a perfect solution to obtain the
  // precision of a field - definition is not clear
  // to be removed
  // ---------

  template <class Field>
  struct Precision;

  template<>
  struct Precision< double >
  {
    static const unsigned int value = 64;
  };

  template<>
  struct Precision< long double >
  {
    static const unsigned int value = 80;
  };

  template<>
  struct Precision< float >
  {
    static const unsigned int value = 32;
  };

#if HAVE_GMP
  template< unsigned int precision >
  struct Precision< GMPField< precision > >
  {
    static const unsigned int value = precision;
  };
#endif

  // ComputeField
  // ------------

  template <class Field,unsigned int sum>
  struct ComputeField
  {
    typedef Field Type;
  };

#if HAVE_GMP
  template< unsigned int precision, unsigned int sum >
  struct ComputeField< GMPField< precision >, sum >
  {
    typedef GMPField<precision+sum> Type;
  };
#endif
} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_UTILITY_FIELD_HH
