// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK_2D_LOCALBASIS_HH
#define DUNE_RANNACHER_TUREK_2D_LOCALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{

  template< class D, class R >
  struct RannacherTurek2DLocalBasis
  {
    typedef LocalBasisTraits< D, 2, FieldVector< D, 2 >,
        R, 1, FieldVector< R, 1 >,
        FieldMatrix< R, 1, 2 > > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 4;
    }

    //! \brief evaluate all shape functions
    inline void evaluateFunction ( const typename Traits::DomainType &in,
                                   std::vector< typename Traits::RangeType > &out ) const
    {
      out.resize(4);
      typename Traits::DomainFieldType qbase = in[0]*in[0]-in[1]*in[1];
      out[0] =  .75 - 2*in[0] +   in[1] + qbase;
      out[1] = -.25           +   in[1] + qbase;
      out[2] =  .75 +   in[0] - 2*in[1] - qbase;
      out[3] = -.25 +   in[0]           - qbase;
    }

    //! \brief evaluate jacobian of all shape functions
    inline void evaluateJacobian ( const typename Traits::DomainType &in,
                                   std::vector< typename Traits::JacobianType > &out ) const
    {
      out.resize(4);

      // see http://www.dune-project.org/doc/doxygen/html/classDune_1_1C1LocalBasisInterface.html#d6f8368f8aa43439cc7ef10419f6e2ea
      // out[i][j][k] = d_k \phi^i_j , where \phi^i_j is the j'th component of the i'th shape function.

      out[0][0][0] = -2 + 2*in[0]; out[0][0][1] =  1 - 2*in[1];
      out[1][0][0] =      2*in[0]; out[1][0][1] =  1 - 2*in[1];
      out[2][0][0] =  1 - 2*in[0]; out[2][0][1] = -2 + 2*in[1];
      out[3][0][0] =  1 - 2*in[0]; out[3][0][1] =      2*in[1];
    }

    //! \brief polynomial order of the shape functions
    unsigned int order () const
    {
      // must be 2 here since it contains x^2 and x^2
      return 2;
    }
  };

} //namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_2D_LOCALBASIS_HH
