// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIINDEX_HH
#define DUNE_MULTIINDEX_HH

#include <dune/common/fvector.hh>

#include <dune/finiteelements/field.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int dim >
  class MultiIndex;

  template< int dim >
  std::ostream &operator<< ( std::ostream &, const MultiIndex< dim > & );



  // MultiIndex
  // ----------

  template< int dim >
  class MultiIndex
  {
    typedef MultiIndex< dim > This;

    friend std::ostream &operator<<<> ( std::ostream &, const This & );

  public:
    static const int dimension = dim;

    MultiIndex ()
      : vecZ_( 0 ),
        vecOMZ_( 0 )
    {}

    MultiIndex ( int, const MultiIndex &other )
      : vecZ_( other.vecOMZ_ ),
        vecOMZ_( other.vecZ_ )
    {}

    MultiIndex ( const This &other )
      : vecZ_( other.vecZ_ ),
        vecOMZ_( other.vecOMZ_ )
    {}

    This &operator= ( const This &other )
    {
      vecZ_   = other.vecZ_;
      vecOMZ_ = other.vecOMZ_;
      return *this;
    }

    This &operator*= ( const This &other )
    {
      vecZ_   += other.vecZ_;
      vecOMZ_ += other.vecOMZ_;
      return *this;
    }

    This &operator/= ( const This &other )
    {
      vecZ_   -= other.vecZ_;
      vecOMZ_ -= other.vecOMZ_;
      return *this;
    }

    This operator* ( const This &other ) const
    {
      This z = *this;
      return (z *= other);
    }

    This operator/ ( const This &other ) const
    {
      This z = *this;
      return (z /= other);
    }

    void set ( int d, int power = 1 )
    {
      vecZ_[ d ] = power;
    }

    int absZ () const
    {
      int ret = 0;
      for( int i = 0; i < dimension; ++i )
        ret += std::abs( vecZ_[ i ] );
      return ret;
    }

    int absOMZ() const
    {
      int ret = 0;
      for( int i = 0; i < dimension; ++i )
        ret += std::abs( vecOMZ_[ i ] );
      return ret;
    }

  private:
    typedef Dune::FieldVector< int, dimension > Vector;

    Vector vecZ_;
    Vector vecOMZ_;
  };


  template< int dim >
  inline std::ostream &
  operator<< ( std::ostream &out, const MultiIndex< dim > &multiIndex )
  {
    if( multiIndex.absZ() == 0 )
      out << "1";
    else
    {
      int absVal = 0;
      for( int i = 0; i < dim; ++i )
      {
        if( multiIndex.vecZ_[ i ] == 0 )
          continue;

        out << char( 'a'+i );
        if( multiIndex.vecZ_[ i ] > 1 )
          out << "**(" << multiIndex.vecZ_[ i ] << ")";
        else if( multiIndex.vecZ_[ i ] < 0 )
          out << "**(" << multiIndex.vecZ_[ i ] << ")";

        absVal += multiIndex.vecZ_[ i ];
        if( absVal < multiIndex.absZ() )
          out << "*";
      }
    }

    if( multiIndex.absOMZ() > 0 )
    {
      for( int i = 0; i < dim; ++i )
      {
        if( multiIndex.vecOMZ_[ i ] == 0 )
          continue;

        out << "(1 - " << char( 'a'+i ) << ")";
        if( multiIndex.vecOMZ_[ i ] > 1 )
          out << "**(" << multiIndex.vecOMZ_[ i ] << ")";
        else if( multiIndex.vecOMZ_[ i ] < 0 )
          out << "**(" << multiIndex.vecOMZ_[ i ] << ")";
      }
    }

    return out;
  }



  template< int dim >
  struct Unity< MultiIndex< dim > >
  {
    typedef MultiIndex< dim > Field;

    operator Field () const
    {
      return Field();
    }

    Field operator- ( const Field &other ) const
    {
      return Field( 1, other );
    }

    Field operator/ ( const Field &other ) const
    {
      return Field() / other;
    }
  };



  template< int dim >
  struct Zero< MultiIndex< dim > >
  {
    typedef MultiIndex< dim > Field;

    // zero does not acutally exist
    operator Field ()
    {
      assert( false );
      return Field();
    }
  };

  template< int dim >
  bool operator< ( const Zero< MultiIndex< dim > > &, const MultiIndex< dim > & )
  {
    return true;
  }

  template< int dim >
  bool operator< ( const MultiIndex< dim > &f, const Zero< MultiIndex< dim > > & )
  {
    return true;
  }

}

#endif // #ifndef DUNE_MULTIINDEX_HH
