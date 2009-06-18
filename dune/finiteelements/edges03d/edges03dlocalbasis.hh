// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALBASIS_HH
#define DUNE_EDGES03DLOCALBASIS_HH

#include <cmath>

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Experimental lowest order edge elements for tetrahedrons.

     (S for simplex)

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class EdgeS03DLocalBasis
    : public C1LocalBasisInterface<
          C1LocalBasisTraits<
              D,3,Dune::FieldVector<D,3>,
              R,3,Dune::FieldVector<R,3>,
              Dune::FieldVector<Dune::FieldVector<R,3>,3>
              >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          , EdgeS03DLocalBasis<D,R>
#endif
          >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<
        D,3,Dune::FieldVector<D,3>,
        R,3,Dune::FieldVector<R,3>,
        Dune::FieldVector<Dune::FieldVector<R,3>,3>
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

    /** \brief Evaluate all shape functions
     *
     * Shape functions are taken from "The finite Element Method in
     * Electromagnetics" by Jianming Jin.
     *
     * In that book the tetrahedron looks like
     * \image html dune/finiteelements/edges03d/Jin2002-reftetrahedron.png
     * Edges have an orientation defined by the order of their nodes:
     * <table>
     * <tr><th>Edge \f$i\f$</th><th>Node \f$i_1\f$</th><th>Node \f$i_2\f$</th></tr>
     * <tr><td>           1</td><td>             1</td><td>             2</td></tr>
     * <tr><td>           2</td><td>             1</td><td>             3</td></tr>
     * <tr><td>           3</td><td>             1</td><td>             4</td></tr>
     * <tr><td>           4</td><td>             2</td><td>             3</td></tr>
     * <tr><td>           5</td><td>             4</td><td>             2</td></tr>
     * <tr><td>           6</td><td>             3</td><td>             4</td></tr>
     * </table>
     * The shape functions are defined as
     * \f[
     *    \mathbf N^e_i=(L^e_{i_1}\nabla L^e_{i_2}
     *                   -L^e_{i_2}\nabla L^e_{i_1})\ell^e_i
     * \f]
     * where \f$\ell^e_i\f$ is the length of edge \f$i\f$, and \f$L^e_j\f$ are
     * the \f$P_1\f$ shape functions:
     * \f[
     *    L^e_j=\frac1{6V^e}(a^e_j+b^e_jx+c^e_jy+d^e_jz)
     * \f]
     * The coefficients \f$a\f$, \f$b\f$, \f$c\f$, and \f$d\f$ can be derived
     * from the following conditions:
     * \f[
     *    L^e_j(x^e_k,y^e_k,z^e_k)=\delta_{jk}\qquad\forall j,k\in\{1\ldots4\}
     * \f]
     * where \f$(x^e_k,y^e_k,z^e_k)\f$ are the coordinates of vertex \f$k\f$.
     *
     * For DUNE we have to replace the indices from Jin with zero-starting
     * indices from DUNE.  We use the following map for the vertex indices:
     * <table>
     * <tr><th>Vertex # (Jin)</th><th>Vertex # (new DUNE)</th></tr>
     * <tr><td>             1</td><td>                  3</td></tr>
     * <tr><td>             2</td><td>                  0</td></tr>
     * <tr><td>             3</td><td>                  1</td></tr>
     * <tr><td>             4</td><td>                  2</td></tr>
     * </table>
     * From that follows the following map for the edge numbers:
     * <table>
     * <tr><th>Edge # (Jin)</th><th>Edge # (new DUNE)</th></tr>
     * <tr><td>           1</td><td>                3</td></tr>
     * <tr><td>           2</td><td>                4</td></tr>
     * <tr><td>           3</td><td>                5</td></tr>
     * <tr><td>           4</td><td>                0</td></tr>
     * <tr><td>           5</td><td>                1</td></tr>
     * <tr><td>           6</td><td>                2</td></tr>
     * </table>
     * For the edge shape functions we will introduce signs
     * \f$s_0,\ldots,s_5\in\{-1,1\}\f$.  Since shape functions of different
     * elements on the same edge must have matching tangential components,
     * these signs must be determined externally in the end.  While indices
     * from Jin have been denoted with \f$i\f$, \f$j\f$ and \f$k\f$ we will
     * denote DUNE indices with lower case greek letters.
     *
     * First we rerwite the mapping from edges to indices into DUNE
     * numbering.  We don't care about the orientation of the edges, since it
     * only influences the sign, which will be determined seperately in the
     * end anyway.
     * <table>
     * <tr><th>Edge \f$\alpha\f$</th>
     *                <th>Node \f$\alpha_1\f$</th>
     *                           <th>Node \f$\alpha_2\f$</th>
     *                                      <th>\f$\ell^e_\alpha\f$</th></tr>
     * <tr><td> 0</td><td> 0</td><td> 1</td><td>\f$1\f$      </td></tr>
     * <tr><td> 1</td><td> 0</td><td> 2</td><td>\f$1\f$      </td></tr>
     * <tr><td> 2</td><td> 1</td><td> 2</td><td>\f$\sqrt2\f$ </td></tr>
     * <tr><td> 3</td><td> 0</td><td> 3</td><td>\f$1\f$      </td></tr>
     * <tr><td> 4</td><td> 1</td><td> 3</td><td>\f$\sqrt2\f$ </td></tr>
     * <tr><td> 5</td><td> 2</td><td> 3</td><td>\f$\sqrt2\f$ </td></tr>
     * </table>
     * The shape functions with DUNE indices are
     * \f[
     *    \mathbf N^e_\alpha=s_\alpha(L^e_{\alpha_1}\nabla L^e_{\alpha_2}
     *                       -L^e_{\alpha_2}\nabla L^e_{\alpha_1})\ell^e_\alpha
     * \f]
     * The corresponding \f$P_1\f$ functions are (with the factor \f$1/6V^e\f$
     * merged into the coefficients).
     * \f[
     *    L^e_\beta=a^e_\beta+b^e_\beta x+c^e_\beta y+d^e_\beta z
     * \f]
     * and the conditions to derive the coefficients are
     * \f[
     *    L^e_\beta(x^e_\gamma,y^e_\gamma,z^e_\gamma)=
     *       \delta_{\beta\gamma}\qquad\forall\beta,\gamma\in\{0\ldots5\}
     * \f]
     * The coordinates of the vertices are
     * \f[ x^e_0=0 \qquad y^e_0=0 \qquad z^e_0=0 \f]
     * \f[ x^e_1=1 \qquad y^e_1=0 \qquad z^e_1=0 \f]
     * \f[ x^e_2=0 \qquad y^e_2=1 \qquad z^e_2=0 \f]
     * \f[ x^e_3=0 \qquad y^e_3=0 \qquad z^e_3=1 \f]
     * Inserting that into the above conditions we arrive at
     * \f[ a^e_0= 1 \qquad b^e_0=-1 \qquad c^e_0=-1 \qquad d^e_0=-1 \f]
     * \f[ a^e_1= 0 \qquad b^e_1= 1 \qquad c^e_1= 0 \qquad d^e_1= 0 \f]
     * \f[ a^e_2= 0 \qquad b^e_2= 0 \qquad c^e_2= 1 \qquad d^e_2= 0 \f]
     * \f[ a^e_3= 0 \qquad b^e_3= 0 \qquad c^e_3= 0 \qquad d^e_3= 1 \f]
     * This yield the following \f$P_1\f$ shape functions and their
     * derivatives:
     * \f[ L^e_0=1-x-y-z \qquad L^e_1=x \qquad L^e_2=y \qquad L^e_3=z \f]
     * \f[
     *    \nabla L^e_0=\begin{pmatrix}-1\\-1\\-1\end{pmatrix} \qquad
     *    \nabla L^e_1=\begin{pmatrix}1\\0\\0\end{pmatrix} \qquad
     *    \nabla L^e_2=\begin{pmatrix}0\\1\\0\end{pmatrix} \qquad
     *    \nabla L^e_3=\begin{pmatrix}0\\0\\1\end{pmatrix}
     * \f]
     * This results in the following edge shape functions:
     * \f[
     *    \mathbf N^e_0=s_0      \begin{pmatrix}1-y-z\\x    \\x    \end{pmatrix} \qquad
     *    \mathbf N^e_1=s_1      \begin{pmatrix}y    \\1-x-z\\y    \end{pmatrix} \qquad
     *    \mathbf N^e_2=s_2\sqrt2\begin{pmatrix}-y   \\x    \\0    \end{pmatrix}
     * \f]
     * \f[
     *    \mathbf N^e_3=s_3      \begin{pmatrix}z    \\z    \\1-x-y\end{pmatrix} \qquad
     *    \mathbf N^e_4=s_4\sqrt2\begin{pmatrix}-z   \\0    \\x    \end{pmatrix} \qquad
     *    \mathbf N^e_5=s_5\sqrt2\begin{pmatrix}0    \\-z   \\y    \end{pmatrix}
     * \f]
     * The Jacobians \f$(J(\mathbf
     * N^e_\alpha))_{\beta\gamma}=\partial_\gamma(\mathbf N^e_\alpha)_\beta\f$
     * look like
     * \f[
     *    J(\mathbf N^e_0)=s_0      \begin{pmatrix}
     *                                  0 & -1 & -1 \\
     *                                  1 &  0 &  0 \\
     *                                  1 &  0 &  0
     *                              \end{pmatrix}   \qquad
     *    J(\mathbf N^e_1)=s_1      \begin{pmatrix}
     *                                  0 &  1 &  0 \\
     *                                 -1 &  0 & -1 \\
     *                                  0 &  1 &  0
     *                              \end{pmatrix}   \qquad
     *    J(\mathbf N^e_2)=s_2\sqrt2\begin{pmatrix}
     *                                  0 & -1 &  0 \\
     *                                  1 &  0 &  0 \\
     *                                  0 &  0 &  0
     *                              \end{pmatrix}
     * \f]
     * \f[
     *    J(\mathbf N^e_3)=s_3      \begin{pmatrix}
     *                                  0 &  0 &  1 \\
     *                                  0 &  0 &  1 \\
     *                                 -1 & -1 &  0
     *                              \end{pmatrix}   \qquad
     *    J(\mathbf N^e_4)=s_4\sqrt2\begin{pmatrix}
     *                                  0 &  0 & -1 \\
     *                                  0 &  0 &  0 \\
     *                                  1 &  0 &  0
     *                              \end{pmatrix}   \qquad
     *    J(\mathbf N^e_5)=s_5\sqrt2\begin{pmatrix}
     *                                  0 &  0 &  0 \\
     *                                  0 &  0 & -1 \\
     *                                  0 &  1 &  0
     *                              \end{pmatrix}
     * \f]
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(6);

      out[0][0] = s[0]*(1      -in[1]-in[2]);
      out[0][1] = s[0]*(  in[0]            );
      out[0][2] = s[0]*(  in[0]            );

      out[1][0] = s[1]*(        in[1]      );
      out[1][1] = s[1]*(1-in[0]      -in[2]);
      out[1][2] = s[1]*(        in[1]      );

      out[2][0] = s[2]*(       -in[1]      );
      out[2][1] = s[2]*(  in[0]            );
      out[2][2] = 0;

      out[3][0] = s[3]*(              in[2]);
      out[3][1] = s[3]*(              in[2]);
      out[3][2] = s[3]*(1-in[0]-in[1]      );

      out[4][0] = s[4]*(             -in[2]);
      out[4][1] = 0;
      out[4][2] = s[4]*(  in[0]            );

      out[5][0] = 0;
      out[5][1] = s[5]*(             -in[2]);
      out[5][2] = s[5]*(        in[1]      );
    }

    //! Get global values of shape functions
    /**
     * Experimental interface for getting global values of shape functions.
     *
     * <b>WARNING: this global interface will result in scaled shape functions
     * compared to the shapefunctions from the local interface with the
     * apropriate transformation applied.  Don't use both interfaces in the
     * same code.</b>
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
     * <b>WARNING: this global interface will result in scaled shape functions
     * compared to the shapefunctions from the local interface with the
     * apropriate transformation applied.  Don't use both interfaces in the
     * same code.</b>
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
      for(int i = 0; i < 3; ++i) {
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

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(6);

      out[0][0][0] =     0; out[0][0][1] = -s[0]; out[0][0][2] = -s[0];
      out[0][1][0] =  s[0]; out[0][1][1] =     0; out[0][1][2] =     0;
      out[0][2][0] =  s[0]; out[0][2][1] =     0; out[0][2][2] =     0;

      out[1][0][0] =     0; out[1][0][1] =  s[1]; out[1][0][2] =     0;
      out[1][1][0] = -s[1]; out[1][1][1] =     0; out[1][1][2] = -s[1];
      out[1][2][0] =     0; out[1][2][1] =  s[1]; out[1][2][2] =     0;

      out[2][0][0] =     0; out[2][0][1] = -s[2]; out[2][0][2] =     0;
      out[2][1][0] =  s[2]; out[2][1][1] =     0; out[2][1][2] =     0;
      out[2][2][0] =     0; out[2][2][1] =     0; out[2][2][2] =     0;

      out[3][0][0] =     0; out[3][0][1] =     0; out[3][0][2] =  s[3];
      out[3][1][0] =     0; out[3][1][1] =     0; out[3][1][2] =  s[3];
      out[3][2][0] = -s[3]; out[3][2][1] = -s[3]; out[3][2][2] =     0;

      out[4][0][0] =     0; out[4][0][1] =     0; out[4][0][2] = -s[4];
      out[4][1][0] =     0; out[4][1][1] =     0; out[4][1][2] =     0;
      out[4][2][0] =  s[4]; out[4][2][1] =     0; out[4][2][2] =     0;

      out[5][0][0] =     0; out[5][0][1] =     0; out[5][0][2] =     0;
      out[5][1][0] =     0; out[5][1][1] =     0; out[5][1][2] = -s[5];
      out[5][2][0] =     0; out[5][2][1] =  s[5]; out[5][2][2] =     0;
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
