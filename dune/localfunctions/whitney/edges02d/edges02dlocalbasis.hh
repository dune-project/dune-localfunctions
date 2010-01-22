// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES02DLOCALBASIS_HH
#define DUNE_EDGES02DLOCALBASIS_HH

#include <cmath>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Experimental lowest order edge elements for triangles.

     (S for simplex)

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

     \note This class does not implement the usual LocalBasisInterface since
           that does not make much sense for vector valued elements.  An
           experimental interface providing global rather than local values is
           provided instead.  Be warned that this interface is subject to
           change without notice, however.

         \nosubgrouping
   */
  template<class D, class R>
  class EdgeS02DLocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<
        D,2,Dune::FieldVector<D,2>,
        R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2>
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

    //! Get global values of shape functions
    /**
     * Experimental interface for getting global values of shape functions.
     *
     * Vector valued shape functions often require a complicated
     * transformation to get the global values from the local values:
     * \f[ \psi_i(g(\mathbf{\hat x}))=L(\hat\psi_i(\mathbf x)) \f]
     * where \f$g\f$ is transformation for coordinates, \f$L\f$ the
     * tranformation for values and the hat \f$\hat{\phantom x}\f$ denotes
     * local quantities.  Worse, if we require for edge elements that
     * \f[ \psi_i\cdot\mathbf t_j=\delta_{ij}\text{ on edge }j \f]
     * and
     * \f[ \hat\psi_i\cdot\mathbf{\hat t}_j=\delta_{ij}\text{ on edge }j \f]
     * then we need a different \f$L\f$ for each shape function.  This
     * property is very desirable be eases debugging and the construction of
     * constraints for grids with hanging nodes.
     *
     * \tparam     Geometry Type of geometry
     * \param[in]  in       Where to evaluate.  <b>NOTE: these are local
     *                      coordinates</b>
     * \param[out] out      The global values of the shape functions.
     * \param[in]  geometry The geometry used for the local to global mapping.
     */
    template<typename Geometry>
    inline void evaluateFunctionGlobal
      (const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out,
      const Geometry &geometry) const
    {
      typename Traits::DomainType pos = geometry.global(in);
      // store coefficients
      FieldVector<typename Traits::DomainFieldType, 3> coeff[3];

      coefficientsGlobal(coeff, geometry);

      out.resize(3);
      for(int i = 0; i < 3; ++i) {
        // \phi^i
        out[i][0] = s[i]*( coeff[i][alpha]*pos[1] + coeff[i][a0]);
        out[i][1] = s[i]*(-coeff[i][alpha]*pos[0] + coeff[i][a1]);
      }
    }

    //! Evaluate global Jacobian of all shape functions
    /**
     * Experimental interface for getting global values of the Jacobian of the
     * shapefunctions.
     *
     * This calculates the global Jacobian
     * \f[
     *   \mathrm J(\psi^i)=
     *   \begin{pmatrix}
     *     \partial\psi^i_0/\partial x_0 & \partial\psi^i_0/\partial x_1 & \cdots \\
     *     \partial\psi^i_1/\partial x_0 & \partial\psi^i_1/\partial x_1 & \cdots \\
     *     \vdots & \vdots & \ddots
     *   \end{pmatrix}
     * \f]
     * Note that this are the derivatives of the global values by the global
     * coordinates, evaluated at local coordinates.
     *
     * \tparam     Geometry Type of geometry
     * \param[in]  in       Where to evaluate.  <b>NOTE: these are local
     *                      coordinates</b>
     * \param[out] out      The values of the global Jacobian.
     * \param[in]  geometry The geometry used for the local to global mapping.
     */
    template<typename Geometry>
    inline void evaluateJacobianGlobal
      (const typename Traits::DomainType& in,         // position
      std::vector<typename Traits::JacobianType>& out,        // return value
      const Geometry &geometry) const
    {
      // store coefficients
      FieldVector<typename Traits::DomainFieldType, 3> coeff[3];

      coefficientsGlobal(coeff, geometry);

      out.resize(3);
      for(int i = 0; i < 3; ++i) {
        // \phi^i
        out[i][0][0] = 0;                     out[i][0][1] = s[i]*coeff[i][alpha];
        out[i][1][0] = -s[i]*coeff[i][alpha]; out[i][1][1] = 0;
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }

  private:
    template<typename Geometry>
    void coefficientsGlobal
      (FieldVector<typename Traits::DomainFieldType, 3> (&coeff)[3],
      const Geometry &geometry) const
    {
      assert(geometry.type().isTriangle());

      /**
       * We have
       * \f[
       *    \psi^i=\alpha^i\mathds{\tilde1}\mathbf x+\mathbf a^i
       * \f]
       * with
       * \f[
       *    \mathds{\tilde1}=\begin{pmatrix}0&1\\-1&0\end{pmatrix}
       * \f]
       * For each shape function we have the equation system
       * \f[
       *    \begin{pmatrix}
       *      (\mathbf x^2-\mathbf x^1)\cdot\mathds{\tilde1}\mathbf x^0
       *              & \mathbf x^2_0-\mathbf x^1_0
       *                        & \mathbf x^2_1-\mathbf x^1_1  \\
       *      (\mathbf x^2-\mathbf x^0)\cdot\mathds{\tilde1}\mathbf x^1
       *              & \mathbf x^2_0-\mathbf x^0_0
       *                        & \mathbf x^2_1-\mathbf x^0_1  \\
       *      (\mathbf x^1-\mathbf x^0)\cdot\mathds{\tilde1}\mathbf x^2
       *              & \mathbf x^1_0-\mathbf x^0_0
       *                        & \mathbf x^1_1-\mathbf x^0_1
       *    \end{pmatrix}
       *    \begin{pmatrix}
       *      \alpha^i      \\
       *      \mathbf a^i_0 \\
       *      \mathbf a^i_1
       *    \end{pmatrix}
       *    =
       *    \begin{pmatrix}
       *      \delta_{i1}\ell^1-\delta_{i0}\ell^0 \\
       *      \delta_{i2}\ell^2+\delta_{i0}\ell^0 \\
       *      \delta_{i1}\ell^1-\delta_{i2}\ell^2
       *    \end{pmatrix}
       * \f]
       * where \f$\mathbf x^i\f$ is the global coordinate of vertex \f$i\f$
       * and \f$\ell^i\f$ is the length (global) of edge \f$i\f$.
       */
      typename Traits::DomainType vertex[3];
      for(int i = 0; i < 3; ++i)
        vertex[i] = geometry.corner(i);

      typename Traits::DomainType offset[3];
      offset[0] = vertex[1]; offset[0] -= vertex[0];
      offset[1] = vertex[2]; offset[1] -= vertex[0];
      offset[2] = vertex[2]; offset[2] -= vertex[1];

      // assemble matrix
      FieldMatrix<typename Traits::DomainFieldType, 3, 3> M;
      M[0][0] = offset[2][0]*vertex[0][1]-offset[2][1]*vertex[0][0];
      M[0][1] = offset[2][0];
      M[0][2] = offset[2][1];

      M[1][0] = offset[1][0]*vertex[1][1]-offset[1][1]*vertex[1][0];
      M[1][1] = offset[1][0];
      M[1][2] = offset[1][1];

      M[2][0] = offset[0][0]*vertex[2][1]-offset[0][1]*vertex[2][0];
      M[2][1] = offset[0][0];
      M[2][2] = offset[0][1];

      M.invert();

      FieldVector<typename Traits::DomainFieldType, 3> rhs;

      // \phi^0
      rhs[0] = -offset[0].two_norm();
      rhs[1] =  offset[0].two_norm();
      rhs[2] = 0;
      M.mv(rhs, coeff[0]);

      // \phi^1
      rhs[0] = offset[1].two_norm();
      rhs[1] = 0;
      rhs[2] = offset[1].two_norm();
      M.mv(rhs, coeff[1]);

      // \phi^2
      rhs[0] = 0;
      rhs[1] =  offset[2].two_norm();
      rhs[2] = -offset[2].two_norm();
      M.mv(rhs, coeff[2]);
    }

    // indices into the coefficient vectors for clearer code
    static const typename FieldVector<typename Traits::DomainFieldType, 3>::size_type
    alpha = 0, a0 = 1, a1 = 2;

    //! The signs
    R s[3];
  };
}
#endif // DUNE_EDGES02DLOCALBASIS_HH
