// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALBASIS_HH
#define DUNE_EDGES03DLOCALBASIS_HH

#include <cmath>

#include <dune/common/fmatrix.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include "../common/localbasis.hh"

namespace Dune
{
  /**\brief Experimental lowest order edge elements for tetrahedrons.

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
  class EdgeS03DLocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<
        D,3,Dune::FieldVector<D,3>,
        R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3>
        > Traits;

    //! contruct a local basis instance with default orientations
    EdgeS03DLocalBasis()
    {
      for(int i = 0; i < 6; ++i)
        s[i] = 1;
    }

    //! contruct a local basis instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS03DLocalBasis(unsigned int orientations)
    {
      for(int i = 0; i < 6; ++i)
        s[i] = 1;
      for(int i = 0; i < 6; ++i)
        if(orientations & (1<<i)) s[i] *= -1;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
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
      FieldVector<typename Traits::DomainFieldType, 6> coeff[6];

      coefficientsGlobal(coeff, geometry);

      out.resize(6);
      for(int i = 0; i < 6; ++i) {
        // \phi^i
        out[i][0] = s[i]*( coeff[i][A01]*pos[1]+coeff[i][A02]*pos[2]+coeff[i][a0]);
        out[i][1] = s[i]*(-coeff[i][A01]*pos[0]+coeff[i][A12]*pos[2]+coeff[i][a1]);
        out[i][2] = s[i]*(-coeff[i][A02]*pos[0]-coeff[i][A12]*pos[1]+coeff[i][a2]);
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
      FieldVector<typename Traits::DomainFieldType, 6> coeff[6];

      coefficientsGlobal(coeff, geometry);

      out.resize(6);
      for(int i = 0; i < 6; ++i) {
        // \phi^i
        out[i][0][0] =  0;
        out[i][0][1] =  s[i]*coeff[i][A01];
        out[i][0][2] =  s[i]*coeff[i][A02];

        out[i][1][0] = -s[i]*coeff[i][A01];
        out[i][1][1] =  0;
        out[i][1][2] =  s[i]*coeff[i][A12];

        out[i][2][0] = -s[i]*coeff[i][A02];
        out[i][2][1] = -s[i]*coeff[i][A12];
        out[i][2][2] =  0;
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
      (FieldVector<typename Traits::DomainFieldType, 6> (&coeff)[6],
      const Geometry &geometry) const
    {
      static const GenericReferenceElement<typename Traits::DomainFieldType, 3> &refElem
        = GenericReferenceElements<typename Traits::DomainFieldType, 3>::simplex();

      assert(geometry.type().isTetrahedron());

      typename Traits::DomainType vertex[4];
      for(int i = 0; i < 4; ++i)
        vertex[i] = geometry.corner(i);

      // assemble matrix and distance vectors
      FieldMatrix<typename Traits::DomainFieldType, 6, 6> M;
      typename Traits::DomainType distance[6];
      for(int j = 0; j < 6; ++j) {
        int v0 = refElem.subEntity(j, 2, 0, 3);
        int v1 = refElem.subEntity(j, 2, 1, 3);
        if(v0 > v1) std::swap(v0, v1);

        distance[j] = vertex[v1]; distance[j] -= vertex[v0];

        M[j][0] = vertex[v1][0]*vertex[v0][1] - vertex[v1][1]*vertex[v0][0];
        M[j][1] = vertex[v1][0]*vertex[v0][2] - vertex[v1][2]*vertex[v0][0];
        M[j][2] = vertex[v1][1]*vertex[v0][2] - vertex[v1][2]*vertex[v0][1];
        M[j][3] = distance[j][0];
        M[j][4] = distance[j][1];
        M[j][5] = distance[j][2];
      }

      M.invert();

      for(int i = 0; i < 6; ++i) {
        FieldVector<typename Traits::DomainFieldType, 6> rhs(0);
        rhs[i] = distance[i].two_norm();

        M.mv(rhs, coeff[i]);
      }
    }

    // indices into the coefficient vectors for clearer code
    enum { A01, A02, A12, a0, a1, a2 };

    R s[6];
  };
}
#endif // DUNE_EDGES03DLOCALBASIS_HH
