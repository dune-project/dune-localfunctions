// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FIELD_HH
#define DUNE_FIELD_HH

#if HAVE_ALGLIB
#include <alglib/amp.h>
#endif
#if HAVE_GMP
#include <dune/finiteelements/generic/math/gmpfield.hh>
#endif

#include <dune/common/fvector.hh>

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
#if HAVE_ALGLIB
  template< unsigned int precision >
  struct Zero< amp::ampf< precision > >
  {
    typedef amp::ampf< precision > Field;
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

#if HAVE_ALGLIB
  template< unsigned int precision >
  inline void field_cast ( const amp::ampf< precision > &f1, double &f2 )
  {
    f2 = f1.toDouble();
  }

  template< unsigned int precision >
  inline void field_cast ( const amp::ampf< precision > &f1, long double &f2 )
  {
    f2 = f1.toDouble();
  }
#endif

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


  template< class F2,class F1 >
  inline F2 field_cast ( const F1 &f1 )
  {
    F2 f2;
    field_cast( f1, f2 );
    return f2;
  }

  template< class F2,class F1,int dim >
  inline Dune::FieldVector<F2,dim> field_cast ( const Dune::FieldVector<F1,dim> &f1 )
  {
    Dune::FieldVector<F2,dim> f2;
    field_cast( f1, f2 );
    return f2;
  }

  template< class F2,class F1,int dim1, int dim2>
  inline Dune::FieldMatrix<F2,dim1,dim2> field_cast ( const Dune::FieldMatrix<F1,dim1,dim2> &f1 )
  {
    Dune::FieldMatrix<F2,dim1,dim2> f2;
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

#if HAVE_ALGLIB
  template< unsigned int precision >
  struct Precision< amp::ampf< precision > >
  {
    static const unsigned int value = precision;
  };
#endif
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

#if HAVE_ALGLIB
  template< unsigned int precision, unsigned int sum >
  struct ComputeField< amp::ampf< precision >, sum >
  {
    typedef amp::ampf<precision+sum> Type;
  };
#endif
#if HAVE_GMP
  template< unsigned int precision, unsigned int sum >
  struct ComputeField< GMPField< precision >, sum >
  {
    typedef GMPField<precision+sum> Type;
  };
#endif
}

// to be moved to different location...
namespace std
{

#if HAVE_ALGLIB
  template< unsigned int precision >
  inline ostream &
  operator<< ( ostream &out,
               const amp::ampf< precision > &value )
  {
    return out << value.toDec();
  }

  template< unsigned int precision >
  inline amp::ampf< precision > sqrt ( const amp::ampf< precision > &a )
  {
    return amp::sqrt( a );
  }

  template< unsigned int precision >
  inline amp::ampf< precision > abs ( const amp::ampf< precision > &a )
  {
    return amp::abs( a );
  }
#endif

#if HAVE_GMP
  template< unsigned int precision >
  inline Dune::GMPField< precision > sqrt ( const Dune::GMPField< precision > &a )
  {
    return ::sqrt( a );
  }

  template< unsigned int precision >
  inline Dune::GMPField< precision > abs ( const Dune::GMPField< precision > &a )
  {
    return ::abs( a );
  }
#endif
}

#endif // #ifndef DUNE_FIELD_HH
