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
      s[0] = 1; s[1] = 1; s[2] = sr2; s[3] = 1; s[4] = sr2; s[5] = sr2;
    }

    //! contruct a local basis instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS03DLocalBasis(unsigned int orientations)
    {
      s[0] = 1; s[1] = 1; s[2] = sr2; s[3] = 1; s[4] = sr2; s[5] = sr2;
      for(int i = 0; i < 3; ++i)
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
    //! square root of 2
    static const R sr2;

    R s[6];
  };
  template<class D, class R>
  const R EdgeS03DLocalBasis<D,R>::sr2 = std::sqrt(R(2.0));
}
#endif // DUNE_EDGES03DLOCALBASIS_HH
