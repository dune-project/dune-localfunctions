// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES02DLOCALBASIS_HH
#define DUNE_EDGES02DLOCALBASIS_HH

#include <cmath>

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Experimental lowest order edge elements for triangles.

     (S for simplex)

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class EdgeS02DLocalBasis
    : public C1LocalBasisInterface<
          C1LocalBasisTraits<
              D,2,Dune::FieldVector<D,2>,
              R,2,Dune::FieldVector<R,2>,
              Dune::FieldVector<Dune::FieldVector<R,2>,2>
              >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          , EdgeS02DLocalBasis<D,R>
#endif
          >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<
        D,2,Dune::FieldVector<D,2>,
        R,2,Dune::FieldVector<R,2>,
        Dune::FieldVector<Dune::FieldVector<R,2>,2>
        > Traits;

    //! contruct a local basis instance with default orientations
    EdgeS02DLocalBasis()
    {
      s[0] = 1; s[1] = 1; s[2] = 1;
    }

    //! contruct a local basis instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS02DLocalBasis(unsigned int orientations)
    {
      s[0] = 1; s[1] = 1; s[2] = 1;
      for(int i = 0; i < 3; ++i)
        if(orientations & (1<<i)) s[i] *= -1;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 3;
    }

    /** \brief Evaluate all shape functions
     *
     * Shape functions are taken from "The finite Element Method in
     * Electromagnetics" by Jianming Jin.
     *
     * In that book the triangle looks like
     * \image html dune/finiteelements/edges02d/Jin2002-reftriangle.png
     *
     * <!--
                3
                /\
               /  \
       edge 3 /    \ edge 2
             /      \
            /        \
          1 ---------- 2
              egde 1
       -->
     * The shape functions are defined as
     *   \f[ \mathbf N^e_1=(L^e_1\nabla L^e_2-L^e_2\nabla L^e_1)\cdot
     *       \ell^e_1 \f]
     * and cycl.  With \f$\ell^e_j\f$ the length of edge \f$j\f$. $L^e_j$ are
     * the shape function for $P_1$ elements:
     *   \f[ L^e_j(x,y)=\frac1{2\Delta^e} (a^e_j+b^e_j\cdot x + c^e_j\cdot y)
     *   \f]
     * where \f$\Delta^e\f$ is the area of the triangle and
     *   \f[ a^e_1=x^e_2 y^e_3-y^e_2 x^e_3 \f]
     *   \f[ b^e_1=y^e_2-y^e_3 \f]
     *   \f[ c^e_1=x^e_3-x^e_2 \f]
     * and cycl, with \f$(x^e_j,y^e_j)\f$ the coordinates of point \f$j\f$.
     *
     * For mapping to DUNEs (new) reference elements, we use the following
     * table
     * <table>
     *   <tr><th>Jin</th><th>DUNE (new)<br>vertices</th><th>DUNE (new)<br>edges</th></tr>
     *   <tr><td> 1 </td><td> 0                    </td><td> 0                 </td></tr>
     *   <tr><td> 2 </td><td> 1                    </td><td> 2                 </td></tr>
     *   <tr><td> 3 </td><td> 2                    </td><td> 1                 </td></tr>
     * </table>
     *
     * Thus, for DUNE indices, we get the following definitions:
     * \f[ \mathbf N^e_0=(L^e_0\nabla L^e_1-L^e_1\nabla L^e_0)\cdot\ell^e_0 \f]
     * \f[ \mathbf N^e_1=(L^e_2\nabla L^e_0-L^e_0\nabla L^e_2)\cdot\ell^e_1 \f]
     * \f[ \mathbf N^e_2=(L^e_1\nabla L^e_2-L^e_2\nabla L^e_1)\cdot\ell^e_2 \f]
     *
     * \f[ L^e_j(x,y)=\frac1{2\Delta^e} (a^e_j+b^e_j\cdot x + c^e_j\cdot y) \f]
     *
     * \f[ a^e_1=x^e_2 y^e_0-y^e_2 x^e_0 \mbox{ and cycl.} \f]
     * \f[ b^e_1=y^e_2-y^e_0 \mbox{ and cycl.} \f]
     * \f[ c^e_1=x^e_0-x^e_2 \mbox{ and cycl.} \f]
     *
     * To break that down onto the reference element, we insert
     * \f[ (x^e_0,y^e_0)=(0,0)\qquad(x^e_1,y^e_1)=(1,0)\qquad(x^e_2,y^e_2)=(0,1) \f]
     * which yields
     * \f[ a^e_0= 1\qquad a^e_1=0\qquad a^e_2=0 \f]
     * \f[ b^e_0=-1\qquad b^e_1=1\qquad b^e_2=0 \f]
     * \f[ c^e_0=-1\qquad c^e_1=0\qquad b^e_2=1 \f]
     * and, together with \f$\Delta^e=1/2\f$
     * \f[ L^e_0(x,y)=1-x-y\qquad \nabla L^e_0=(-1,-1) \f]
     * \f[ L^e_1(x,y)=x    \qquad \nabla L^e_1=(1,0) \f]
     * \f[ L^e_2(x,y)=y    \qquad \nabla L^e_2=(0,1) \f]
     * Thus we arrive at
     * With \f$\ell^e_0=1\f$, \f$\ell^e_1=1\f$ and \f$\ell^e_2=\sqrt2\f$ we
     * arrive at
     * \f[ \mathbf N^e_0(x,y)=\ell^e_0\cdot(1-y,   x) \f]
     * \f[ \mathbf N^e_1(x,y)=\ell^e_1\cdot( -y,-1+x) \f]
     * \f[ \mathbf N^e_2(x,y)=\ell^e_2\cdot( -y,   x) \f]
     *
     * It turns out that having the factor \f$\ell^e_i\f$ make the
     * transformation of the values from reference to world coordinates much
     * more complicated (if at all possible by a linear transformation).
     * Therefore, we leave this factor out.
     * \f[
     *    \mathbf N^e_i\cdot\mathbf{\hat t}_j=\delta_{ij}\text{ on }\Gamma_j
     * \f]
     * then becomes
     * \f[
     *    \mathbf N^e_i\cdot\mathbf{\hat t}_j=\delta_{ij}/\ell^e_j\text{ on }\Gamma_j
     * \f]
     *
     * These three base functions all turn counterclockwise the viewed in
     * paraview (x-axis pointing right, y-axis pointing up).  When two
     * triangles are joined together, this will lead to the tangential
     * components of the two corresponding shape functions from the two
     * triangles having an opposite sign.  To prevent that, we do two things:
     *  -# we multiply \f$\mathbf N^e_1\f$ by \f$-1\f$ so that every shape
     *     function in the reference element per default points from the
     *     vertexwith the lower index to the vertex with the higher index on
     *     it's corresponding edge, and
     *  -# for each shape function \f$\mathbf N^e_alpha\f$ we introduce a sign
     *     \f$s_\alpha\f$ which can be set externally, e.g. by the
     *     FiniteElementMap from dune-pdelab.
     *  .
     * Thus the final shape functions look like this:
     * \f[ \mathbf N^e_0(x,y)=s_0(1-y,  x) \f]
     * \f[ \mathbf N^e_1(x,y)=s_1(  y,1-x) \f]
     * \f[ \mathbf N^e_2(x,y)=s_2( -y,  x) \f]
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(3);
      out[0][0] = s[0]*(1-in[1]); out[0][1] = s[0]*(  in[0]);
      out[1][0] = s[1]*(  in[1]); out[1][1] = s[1]*(1-in[0]);
      out[2][0] = s[2]*( -in[1]); out[2][1] = s[2]*(  in[0]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(3);
      // basis function 0
      out[0][0][0] =     0; out[0][1][0] = s[0];
      out[0][0][1] = -s[0]; out[0][1][1] =    0;
      // basis function 1
      out[1][0][0] =     0; out[1][1][0] = s[1];
      out[1][0][1] = -s[1]; out[1][1][1] =    0;
      // basis function 2
      out[2][0][0] =     0; out[2][1][0] = s[2];
      out[2][0][1] = -s[2]; out[2][1][1] =    0;
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }

  private:
    //! The signs
    R s[3];
  };
}
#endif // DUNE_EDGES02DLOCALBASIS_HH
