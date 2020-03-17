// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALCOMPUTE_HH
#define DUNE_ORTHONORMALCOMPUTE_HH

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/lfematrix.hh>
#include <dune/localfunctions/utility/monomialbasis.hh>
#include <dune/localfunctions/utility/multiindex.hh>

namespace ONBCompute
{

  template< class scalar_t >
  scalar_t factorial( int start, int end )
  {
    scalar_t ret( 1 );
    for( int j = start; j <= end; ++j )
      ret *= scalar_t( j );
    return ret;
  }



  // Integral
  // --------

  template< class Topology >
  struct Integral;

  template< class Base >
  struct Integral< Dune::Impl::Pyramid< Base > >
  {
    template< int dim, class scalar_t >
    static int compute ( const Dune::MultiIndex< dim, scalar_t > &alpha,
                         scalar_t &p, scalar_t &q )
    {
      const int dimension = Base::dimension+1;
      int i = alpha.z( Base::dimension );
      int ord = Integral< Base >::compute( alpha, p, q );
      p *= factorial< scalar_t >( 1, i );
      q *= factorial< scalar_t >( dimension + ord, dimension + ord + i );
      return ord + i;
    }
  };

  template< class Base >
  struct Integral< Dune::Impl::Prism< Base > >
  {
    template< int dim, class scalar_t >
    static int compute ( const Dune::MultiIndex< dim, scalar_t > &alpha,
                         scalar_t &p, scalar_t &q )
    {
      int i = alpha.z( Base::dimension );
      int ord = Integral< Base >::compute( alpha, p, q );
      //Integral< Base >::compute( alpha, p, q );
      //p *= scalar_t( 1 );
      q *= scalar_t( i+1 );
      return ord + i;
    }
  };

  template<>
  struct Integral< Dune::Impl::Point >
  {
    template< int dim, class scalar_t >
    static int compute ( const Dune::MultiIndex< dim, scalar_t > &alpha,
                         scalar_t &p, scalar_t &q )
    {
      p = scalar_t( 1 );
      q = scalar_t( 1 );
      return 0;
    }
  };



  // ONBMatrix
  // ---------

  template< class Topology, class scalar_t >
  class ONBMatrix
    : public Dune::LFEMatrix< scalar_t >
  {
    typedef ONBMatrix< Topology, scalar_t > This;
    typedef Dune::LFEMatrix< scalar_t > Base;

  public:
    typedef std::vector< scalar_t > vec_t;
    typedef Dune::LFEMatrix< scalar_t > mat_t;

    explicit ONBMatrix ( unsigned int order )
    {
      // get all multiindecies for monomial basis
      const unsigned int dim = Topology::dimension;
      typedef Dune::MultiIndex< dim, scalar_t > MI;
      Dune::StandardMonomialBasis< dim, MI > basis( order );
      const std::size_t size = basis.size();
      std::vector< Dune::FieldVector< MI, 1 > > y( size );
      Dune::FieldVector< MI, dim > x;
      for( unsigned int i = 0; i < dim; ++i )
        x[ i ].set( i );
      basis.evaluate( x, y );

      // set bounds of data
      Base::resize( size, size );
      S.resize( size, size );
      d.resize( size );

      // setup matrix for bilinear form x^T S y: S_ij = int_A x^(i+j)
      scalar_t p, q;
      for( std::size_t i = 0; i < size; ++i )
      {
        for( std::size_t j = 0; j < size; ++j )
        {
          Integral< Topology >::compute( y[ i ][ 0 ] * y[ j ][ 0 ], p, q );
          S( i, j ) = p;
          S( i, j ) /= q;
        }
      }

      // orthonormalize
      gramSchmidt();
    }

    template< class Vector >
    void row ( unsigned int row, Vector &vec ) const
    {
      // transposed matrix is required
      assert( row < Base::cols() );
      for( std::size_t i = 0; i < Base::rows(); ++i )
        Dune::field_cast( Base::operator()( i, row ), vec[ i ] );
    }

  private:
    void sprod ( int col1, int col2, scalar_t &ret )
    {
      ret = 0;
      for( int k = 0; k <= col1; ++k )
      {
        for( int l = 0; l <=col2; ++l )
          ret += Base::operator()( l, col2 ) * S( l, k ) * Base::operator()( k, col1 );
      }
    }

    void vmul ( std::size_t col, std::size_t rowEnd, const scalar_t &s )
    {
      for( std::size_t i = 0; i <= rowEnd; ++i )
        Base::operator()( i, col ) *= s;
    }

    void vsub ( std::size_t coldest, std::size_t colsrc, std::size_t rowEnd, const scalar_t &s )
    {
      for( std::size_t i = 0; i <= rowEnd; ++i )
        Base::operator()( i, coldest ) -= s * Base::operator()( i, colsrc );
    }

    void gramSchmidt ()
    {
      using std::sqrt;
      // setup identity
      const std::size_t N = Base::rows();
      for( std::size_t i = 0; i < N; ++i )
      {
        for( std::size_t j = 0; j < N; ++j )
          Base::operator()( i, j ) = scalar_t( i == j ? 1 : 0 );
      }

      // perform Gram-Schmidt procedure
      scalar_t s;
      sprod( 0, 0, s );
      vmul( 0, 0, scalar_t( 1 ) / sqrt( s ) );
      for( std::size_t i = 1; i < N; ++i )
      {
        for( std::size_t k = 0; k < i; ++k )
        {
          sprod( i, k, s );
          vsub( i, k, i, s );
        }
        sprod( i, i, s );
        vmul( i, i, scalar_t( 1 ) / sqrt( s ) );
      }
    }

    vec_t d;
    mat_t S;
  };

} // namespace ONBCompute

#endif // #ifndef DUNE_ORTHONORMALCOMPUTE_HH
