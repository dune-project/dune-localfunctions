// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES02DLOCALINTERPOLATION_HH
#define DUNE_EDGES02DLOCALINTERPOLATION_HH

#include <cmath>

#include "../common/localinterpolation.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Interpolation for experimental lowest order edge elements for triangles.

         \tparam LB LocalBasisImplementation
         \nosubgrouping
   */

  template<class LB>
  class EdgeS02DLocalInterpolation
    : public LocalInterpolationInterface<EdgeS02DLocalInterpolation<LB> >
  {
  public:
    //! contruct an interpolation instance with default orientations
    EdgeS02DLocalInterpolation()
    {
      s[0] = s[1] = s[2] = 1;
    }

    //! contruct an interpolation instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS02DLocalInterpolation(unsigned int orientations)
    {
      s[0] = s[1] = s[2] = 1;
      for(int i = 0; i < 3; ++i)
        if(orientations & (1<<i)) s[i] = -1;
    }

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(3);

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
       *     \mathbf{\hat t}_0=(1,0)\qquad
       *     \mathbf{\hat t}_1=(0,1)\qquad
       *     \mathbf{\hat t}_2=\frac1{\sqrt2}\cdot(-1,1)
       *  \f]
       *  The scalar product \f$\mathbf N^e_\alpha\cdot\mathbf{\hat
       *  t}_\alpha\f$ then becomes:
       *  \f[
       *     \mathbf N^e_0\cdot\mathbf{\hat t}_0=s_0( 1  -y)\stackrel{  y=0}= s_0 \qquad
       *     \mathbf N^e_1\cdot\mathbf{\hat t}_1=s_1(-1+x  )\stackrel{x  =0}=-s_1 \qquad
       *     \mathbf N^e_2\cdot\mathbf{\hat t}_2=s_2(   x+y)\stackrel{x+y=1}= s_2
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
       *  tangential components on the edges, and were integrating the
       *  tangential components on the edges, but the Gaussian quadrature rule
       *  for constant functions is the same as for first order
       *  polynomials...).  So we can replace the integral
       *  \f$\int_{C_\alpha}ds\,\mathbf a\f$ by \f$\ell^e_\alpha\mathbf
       *  a(x^e_\alpha,y^e_\alpha)\f$:
       *  \f[
       *     (\mathbf a,\mathbf N^e_\alpha)=\ell^e_\alpha
       *     \overline{(\mathbf N^e_\alpha\cdot\mathbf{\hat t}_\alpha)}
       *     (\mathbf a(x^e_\alpha,y^e_\alpha)\cdot\mathbf{\hat t}_\alpha)
       *  \f]
       *  Written out for the different base functions this becomes:
       *  \f[
       *     (\mathbf a,\mathbf N^e_0)=
       *             s_0\mathbf a_x(0.5,0)     \qquad
       *     (\mathbf a,\mathbf N^e_1)=
       *            -s_1\mathbf a_y(0,0.5)     \qquad
       *     (\mathbf a,\mathbf N^e_2)=
       *             s_2(-\mathbf a_x(0.5,0.5)+\mathbf b_y(0.5,0.5))
       *  \f]
       *
       *  Now about the mass matrix \f$M_{\alpha,\beta}=(\mathbf
       *  N^e_\alpha,\mathbf N^e_\beta)\f$.  Since the tangential component of
       *  all base functions are 0 on all edges but one, \f$M\f$ is diagonal.
       *  The components on the diagonal are:
       *  \f[
       *     M_{00}=1 \qquad M_{11}=1 \qquad M_{22}=\sqrt2
       *  \f]
       *  The corresponding coefficients of \f$M^{-1}\f$ are then
       *  \f[
       *     (M^{-1})_{00}=1 \qquad (M^{-1})_{11}=1 \qquad (M^{-1})_{22}=\frac1{\sqrt2}
       *  \f]
       *  Thus we arrive at
       *  \f[
       *     a_0= s_0\mathbf a_x(0.5,0)     \qquad
       *     a_1=-s_1\mathbf a_y(0,0.5)     \qquad
       *     a_2=\frac{s_2}{\sqrt2}
       *        (-\mathbf a_x(0.5,0.5)+\mathbf a_y(0.5,0.5))
       *  \f]
       */

      x[0] = 0.5; x[1] = 0.0; f.evaluate(x,y); out[0] = s[0]*        y[0];
      x[0] = 0.0; x[1] = 0.5; f.evaluate(x,y); out[1] = s[1]*            -y[1];
      x[0] = 0.5; x[1] = 0.5; f.evaluate(x,y); out[2] = s[2]*sr0_5*(-y[0]+y[1]);
    }

  private:
    //! square root of 1/2
    static const typename LB::Traits::RangeFieldType sr0_5;

    typename LB::Traits::RangeFieldType s[3];
  };

  template<class LB>
  const typename LB::Traits::RangeFieldType EdgeS02DLocalInterpolation <LB>::sr0_5 =
    std::sqrt(typename LB::Traits::RangeFieldType(0.5));
}

#endif // DUNE_EDGES02DLOCALINTERPOLATION_HH
