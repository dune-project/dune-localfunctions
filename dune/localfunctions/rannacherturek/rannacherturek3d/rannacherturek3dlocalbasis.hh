// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RANNACHER_TUREK_3D_LOCALBASIS_HH
#define DUNE_RANNACHER_TUREK_3D_LOCALBASIS_HH

#include <numeric>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{

  template< class D, class R >
  class RannacherTurek3DLocalBasis
  {
    static const int coefficients[ 6 ][ 6 ];

  public:
    typedef LocalBasisTraits< D, 3, FieldVector< D, 3 >,
        R, 1, FieldVector< R, 1 >,
        FieldMatrix< R, 1, 3 > > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
    }

    //! \brief evaluate all shape functions
    inline void evaluateFunction ( const typename Traits::DomainType &in,
                                   std::vector< typename Traits::RangeType > &out ) const
    {
      typedef typename Traits::RangeFieldType RangeFieldType;
      RangeFieldType y[ 6 ] = { 1, in[ 0 ], in[ 1 ], in[ 2 ],
                                in[ 0 ]*in[ 0 ] - in[ 1 ]*in[ 1 ],
                                in[ 1 ]*in[ 1 ] - in[ 2 ]*in[ 2 ] };
      out.resize( size() );
      for( unsigned int i = 0; i < size(); ++i )
      {
        out[ i ] = RangeFieldType( 0 );
        for( unsigned int j = 0; j < 6; ++j )
          out[ i ] += coefficients[ i ][ j ]*y[ j ];
        out[ i ] /= RangeFieldType( 3 );
      }
    }

    //! \brief evaluate jacobian of all shape functions
    inline void evaluateJacobian ( const typename Traits::DomainType &in,
                                   std::vector< typename Traits::JacobianType > &out ) const
    {
      typedef typename Traits::RangeFieldType RangeFieldType;
      RangeFieldType y0[ 5 ] = { 1, 0, 0, 2*in[ 0 ], 0 };
      RangeFieldType y1[ 5 ] = { 0, 1, 0, -2*in[ 1 ], 2*in[ 1 ] };
      RangeFieldType y2[ 5 ] = { 0, 0, 1, 0, -2*in[ 2 ] };

      out.resize( size() );
      for( unsigned int i = 0; i < size(); ++i )
      {
        out[ i ] = RangeFieldType( 0 );
        for( unsigned int j = 0; j < 5; ++j )
        {
          out[ i ][ 0 ][ 0 ] += coefficients[ i ][ j+1 ]*y0[ j ];
          out[ i ][ 0 ][ 1 ] += coefficients[ i ][ j+1 ]*y1[ j ];
          out[ i ][ 0 ][ 2 ] += coefficients[ i ][ j+1 ]*y2[ j ];
        }
        out[ i ] /= RangeFieldType( 3 );
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        using RangeFieldType = typename Traits::RangeFieldType;
        RangeFieldType y[3][5] = { { 1.0, 0.0, 0.0,  2*in[0],      0.0 },
                                   { 0.0, 1.0, 0.0, -2*in[1],  2*in[1] },
                                   { 0.0, 0.0, 1.0,      0.0, -2*in[2] } };

        for (std::size_t i = 0; i < size(); ++i) {
          out[i] = RangeFieldType{0};
          for (std::size_t j = 0; j < 5; ++j)
            out[i] += coefficients[i][j+1] * y[direction][j];
          out[i] /= RangeFieldType{3};
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }
  };



  // RannacherTurek3DLocalBasis::coefficients
  // ----------------------------------------

  template< class D, class R >
  const int RannacherTurek3DLocalBasis< D, R >
  ::coefficients[ 6 ][ 6 ] = {{  2, -7,  2,  2,  4,  2 },
                              { -1, -1,  2,  2,  4,  2 },
                              {  2,  2, -7,  2, -2,  2 },
                              { -1,  2, -1,  2, -2,  2 },
                              {  2,  2,  2, -7, -2, -4 },
                              { -1,  2,  2, -1, -2, -4 }};

} //namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_3D_LOCALBASIS_HH
