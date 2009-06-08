// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALINTERPOLATION_HH
#define DUNE_EDGES03DLOCALINTERPOLATION_HH

#include <cmath>

#include "../common/localinterpolation.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Interpolation for experimental lowest order edge elements for tetrahedrons.

         \tparam LB LocalBasisImplementation
         \nosubgrouping
   */

  template<class LB>
  class EdgeS03DLocalInterpolation
    : public LocalInterpolationInterface<EdgeS03DLocalInterpolation<LB> >
  {
  public:
    //! contruct an interpolation instance with default orientations
    EdgeS03DLocalInterpolation()
    {
      s[0] = 1; s[1] = 1; s[2] = sr0_5; s[3] = 1; s[4] = sr0_5; s[5] = sr0_5;
    }

    //! contruct an interpolation instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS03DLocalInterpolation(unsigned int orientations)
    {
      s[0] = 1; s[1] = 1; s[2] = sr0_5; s[3] = 1; s[4] = sr0_5; s[5] = sr0_5;
      for(int i = 0; i < 3; ++i)
        if(orientations & (1<<i)) s[i] *= -1;
    }

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(6);

      /**
       *  A function \f$\mathbf a\f$ can be expressed in terms of a basis
       *  \f$\mathbf N^e_\alpha\f$ as
       *  \f[
       *     \mathbf a=\sum_\alpha a_\alpha\mathbf N^e_\alpha
       *  \f]
       *  with the coefficients \f$a_\alpha\f$, which can be determined as
       *  \f[
       *     a_\alpha=\sum_\beta(\mathbf a,\mathbf N^e_\beta)(M^{-1})_{\beta\alpha}
       *  \f]
       *
       *  The scalar product used for interpolating this element is
       *  \f[
       *     (\mathbf a,\mathbf b)=\sum_\alpha\int_{C_\alpha}ds\,
       *     (\mathbf a\cdot\mathbf{\hat t}_\alpha)
       *     \overline{(\mathbf b\cdot\mathbf{\hat t}_\alpha)}
       *  \f]
       *  where \f$C_\alpha\f$ denotes edge \f$\alpha\f$ and \f$\mathbf{\hat
       *  t}_\alpha\f$ a unitvector tangential to \f$C_\alpha\f$.  The
       *  orientation of \f$\mathbf{\hat t}_\alpha\f$ does not matter, since a
       *  factor of -1 introduced in the first factor under the itegral will
       *  be canceled out by the second.  Here we use the following definition
       *  for \f$\mathbf{\hat t}_\alpha\f$:
       *  \f[
       *     \mathbf{\hat t}_0=              \begin{pmatrix} 1\\ 0\\ 0\end{pmatrix}\qquad
       *     \mathbf{\hat t}_1=              \begin{pmatrix} 0\\ 1\\ 0\end{pmatrix}\qquad
       *     \mathbf{\hat t}_2=\frac1{\sqrt2}\begin{pmatrix}-1\\ 1\\ 0\end{pmatrix}\qquad
       *     \mathbf{\hat t}_3=              \begin{pmatrix} 0\\ 0\\ 1\end{pmatrix}\qquad
       *     \mathbf{\hat t}_4=\frac1{\sqrt2}\begin{pmatrix}-1\\ 0\\ 1\end{pmatrix}\qquad
       *     \mathbf{\hat t}_5=\frac1{\sqrt2}\begin{pmatrix} 0\\-1\\ 1\end{pmatrix}\qquad
       *  \f]
       *  The scalar product \f$\mathbf N^e_\alpha\cdot\mathbf{\hat
       *  t}_\alpha\f$ then becomes:
       *  \f[
       *     \mathbf N^e_0\cdot\mathbf{\hat t}_0=s_0( 1  -y-z)
       *                    \stackrel{ y=0 }{\stackrel{z=0}=}s_0 \qquad
       *     \mathbf N^e_1\cdot\mathbf{\hat t}_1=s_1( 1-x  -z)
       *                    \stackrel{ x=0 }{\stackrel{z=0}=}s_1 \qquad
       *     \mathbf N^e_2\cdot\mathbf{\hat t}_2=s_2(   x+y  )
       *                    \stackrel{x+y=1}{\stackrel{z=0}=}s_2
       *  \f]
       *  \f[
       *     \mathbf N^e_3\cdot\mathbf{\hat t}_3=s_3( 1-x-y  )
       *                    \stackrel{ x=0 }{\stackrel{y=0}=}s_3 \qquad
       *     \mathbf N^e_4\cdot\mathbf{\hat t}_4=s_4(   x  +z)
       *                    \stackrel{x+z=1}{\stackrel{y=0}=}s_4 \qquad
       *     \mathbf N^e_5\cdot\mathbf{\hat t}_5=s_5(     y+z)
       *                    \stackrel{y+z=1}{\stackrel{x=0}=}s_5
       *  \f]
       *  Here the relation which define \f$C_\alpha\f$ have been used in the
       *  second step.  The second factor can then be taken out ouf the
       *  integral.  Further, on all edges \f$C_\beta\f$,
       *  \f$\beta\not=\alpha\f$ \f$\mathbf N^e_\alpha\cdot\mathbf{\hat
       *  t}_\alpha=0\f$.  Thus the sum is comprised only of one summand:
       *  \f[
       *     (\mathbf a,\mathbf N^e_\alpha)=
       *     \overline{(\mathbf N^e_\alpha\cdot\mathbf{\hat t}_\alpha)}
       *     \int_{C_\alpha}ds\,(\mathbf a\cdot\mathbf{\hat t}_\alpha)
       *  \f]
       *  Integration can be exchanged with the second scalar product:
       *  \f[
       *     (\mathbf a,\mathbf N^e_\alpha)=
       *     \overline{(\mathbf N^e_\alpha\cdot\mathbf{\hat t}_\alpha)}
       *     \overline{\left(\mathbf{\hat t}_\alpha\cdot
       *                     \int_{C_\alpha}ds\,\mathbf a\right)}
       *  \f]
       *
       *  For integration we use a hardcoded Gaussian quadrature rule.  Since
       *  we want to interpolate into a linear basis, we use the exact rule
       *  for first order polynomials (actually, the basis is constant in it's
       *  tangential components on the edges, and we're integrating the
       *  tangential components on the edges, but the Gaussian quadrature rule
       *  for constant functions is the same as for first order
       *  polynomials...).  So we can replace the integral
       *  \f$\int_{C_\alpha}ds\,\mathbf a\f$ by \f$\ell^e_\alpha\mathbf
       *  a(x^e_\alpha,y^e_\alpha,z^e_\alpha)\f$:
       *  \f[
       *     (\mathbf a,\mathbf N^e_\alpha)=\ell^e_\alpha
       *     \overline{(\mathbf N^e_\alpha\cdot\mathbf{\hat t}_\alpha)}
       *     \left(\mathbf a\begin{pmatrix}x^e_\alpha\\y^e_\alpha\\z^e_\alpha
       *     \end{pmatrix}\cdot\mathbf{\hat t}_\alpha\right)
       *  \f]
       *  Written out for the different base functions this becomes:
       *  \f[
       *     (\mathbf a,\mathbf N^e_0)=
       *       s_0\mathbf a_x\begin{pmatrix}0.5\\0\\0\end{pmatrix} \qquad
       *     (\mathbf a,\mathbf N^e_1)=
       *       s_1\mathbf a_y\begin{pmatrix}0\\0.5\\0\end{pmatrix} \qquad
       *     (\mathbf a,\mathbf N^e_2)=
       *       s_2\left\{-\mathbf a_x\begin{pmatrix}0.5\\0.5\\0\end{pmatrix}
       *                 +\mathbf a_y\begin{pmatrix}0.5\\0.5\\0\end{pmatrix}
       *          \right\}
       *  \f]
       *  \f[
       *     (\mathbf a,\mathbf N^e_3)=
       *       s_3\mathbf a_z\begin{pmatrix}0\\0\\0.5\end{pmatrix} \qquad
       *     (\mathbf a,\mathbf N^e_4)=
       *       s_4\left\{-\mathbf a_x\begin{pmatrix}0.5\\0\\0.5\end{pmatrix}
       *                 +\mathbf a_z\begin{pmatrix}0.5\\0\\0.5\end{pmatrix}
       *          \right\}                                         \qquad
       *     (\mathbf a,\mathbf N^e_5)=
       *       s_5\left\{-\mathbf a_y\begin{pmatrix}0\\0.5\\0.5\end{pmatrix}
       *                 +\mathbf a_z\begin{pmatrix}0\\0.5\\0.5\end{pmatrix}
       *          \right\}
       *  \f]
       *
       *  Now about the mass matrix \f$M_{\alpha,\beta}=(\mathbf
       *  N^e_\alpha,\mathbf N^e_\beta)\f$.  Since the tangential component of
       *  all base functions are 0 on all edges but one, \f$M\f$ is diagonal.
       *  The components on the diagonal are:
       *  \f[
       *     M=\begin{pmatrix}
       *         1 & 0 &    0   & 0 &    0   &    0   \\
       *         0 & 1 &    0   & 0 &    0   &    0   \\
       *         0 & 0 & \sqrt2 & 0 &    0   &    0   \\
       *         0 & 0 &    0   & 1 &    0   &    0   \\
       *         0 & 0 &    0   & 0 & \sqrt2 &    0   \\
       *         0 & 0 &    0   & 0 &    0   & \sqrt2
       *       \end{pmatrix}
       *  \f]
       *  The corresponding coefficients of \f$M^{-1}\f$ are then
       *  \f[
       *     M^{-1}=\begin{pmatrix}
       *         1 & 0 &     0    & 0 &     0    &     0    \\
       *         0 & 1 &     0    & 0 &     0    &     0    \\
       *         0 & 0 & 1/\sqrt2 & 0 &     0    &     0    \\
       *         0 & 0 &     0    & 1 &     0    &     0    \\
       *         0 & 0 &     0    & 0 & 1/\sqrt2 &     0    \\
       *         0 & 0 &     0    & 0 &     0    & 1/\sqrt2
       *       \end{pmatrix}
       *  \f]
       *  Thus we arrive at
       *  \f[
       *     a_0=s_0
       *       \mathbf a_x\begin{pmatrix}0.5\\0\\0\end{pmatrix} \qquad
       *     a_1=s_1
       *       \mathbf a_y\begin{pmatrix}0\\0.5\\0\end{pmatrix} \qquad
       *     a_2=\frac{s_2}{\sqrt2}
       *       \left\{-\mathbf a_x\begin{pmatrix}0.5\\0.5\\0\end{pmatrix}
       *              +\mathbf a_y\begin{pmatrix}0.5\\0.5\\0\end{pmatrix}
       *       \right\}
       *  \f]
       *  \f[
       *     a_3=s_3
       *       \mathbf a_z\begin{pmatrix}0\\0\\0.5\end{pmatrix} \qquad
       *     a_4=\frac{s_4}{\sqrt2}
       *       \left\{-\mathbf a_x\begin{pmatrix}0.5\\0\\0.5\end{pmatrix}
       *              +\mathbf a_z\begin{pmatrix}0.5\\0\\0.5\end{pmatrix}
       *       \right\}                                         \qquad
       *     a_5=\frac{s_5}{\sqrt2}
       *       \left\{-\mathbf a_y\begin{pmatrix}0\\0.5\\0.5\end{pmatrix}
       *              +\mathbf a_z\begin{pmatrix}0\\0.5\\0.5\end{pmatrix}
       *       \right\}
       *  \f]
       */

      x[0] = 0.5; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y);
      out[0] = s[0]*( y[0]          );

      x[0] = 0.0; x[1] = 0.5; x[2] = 0.0; f.evaluate(x,y);
      out[1] = s[1]*(      y[1]     );

      x[0] = 0.5; x[1] = 0.5; x[2] = 0.0; f.evaluate(x,y);
      out[2] = s[2]*(-y[0]+y[1]     );

      x[0] = 0.0; x[1] = 0.0; x[2] = 0.5; f.evaluate(x,y);
      out[3] = s[3]*(           y[2]);

      x[0] = 0.5; x[1] = 0.0; x[2] = 0.5; f.evaluate(x,y);
      out[4] = s[4]*(-y[0]     +y[2]);

      x[0] = 0.0; x[1] = 0.5; x[2] = 0.5; f.evaluate(x,y);
      out[5] = s[5]*(     -y[1]+y[2]);
    }

  private:
    //! square root of 1/2
    static const typename LB::Traits::RangeFieldType sr0_5;

    typename LB::Traits::RangeFieldType s[6];
  };

  template<class LB>
  const typename LB::Traits::RangeFieldType EdgeS03DLocalInterpolation <LB>::sr0_5 =
    std::sqrt(typename LB::Traits::RangeFieldType(0.5));
}

#endif // DUNE_EDGES03DLOCALINTERPOLATION_HH
