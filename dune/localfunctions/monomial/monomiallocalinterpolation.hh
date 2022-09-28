// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  template<class LB, unsigned int size>
  class MonomialLocalInterpolation
  {
    typedef typename LB::Traits::DomainType D;
    typedef typename LB::Traits::DomainFieldType DF;
    static const int dimD=LB::Traits::dimDomain;
    typedef typename LB::Traits::RangeType R;
    typedef typename LB::Traits::RangeFieldType RF;

    typedef QuadratureRule<DF,dimD> QR;
    typedef typename QR::iterator QRiterator;

  public:
    MonomialLocalInterpolation (const GeometryType &gt_,
                             const LB &lb_)
      : gt(gt_), lb(lb_), Minv(0)
        , qr(QuadratureRules<DF,dimD>::rule(gt, 2*lb.order()))
    {
      // Compute inverse of the mass matrix of the local basis, and store it in Minv
      if(size != lb.size())
        DUNE_THROW(Exception, "size template parameter does not match size of "
                   "local basis");

      const QRiterator qrend = qr.end();
      for(QRiterator qrit = qr.begin(); qrit != qrend; ++qrit) {
        std::vector<R> base;
        lb.evaluateFunction(qrit->position(),base);

        for(unsigned int i = 0; i < size; ++i)
          for(unsigned int j = 0; j < size; ++j)
            Minv[i][j] += qrit->weight() * base[i] * base[j];
      }
      Minv.invert();
    }

    /** \brief Determine coefficients interpolating a given function
     *
     * The method computes the coefficients
     * for the L^2 projection with respect to the given
     * GeometryType. Be careful: the implementation is
     * unstable for higher polynomial degrees.
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      using DomainType = std::decay_t<decltype(qr.begin()->position())>;

      auto&& f = Impl::makeFunctionWithCallOperator<DomainType>(ff);

      out.clear();
      out.resize(size, 0);

      const QRiterator qrend = qr.end();
      for(QRiterator qrit = qr.begin(); qrit != qrend; ++qrit) {
        //TODO: mass matrix
        R y = f(qrit->position());

        std::vector<R> base;
        lb.evaluateFunction(qrit->position(),base);

        for(unsigned int i = 0; i < size; ++i)
          for(unsigned int j = 0; j < size; ++j)
            out[i] += Minv[i][j] * qrit->weight() * y * base[j];
      }
    }

  private:
    GeometryType gt;
    const LB &lb;
    FieldMatrix<RF, size, size> Minv;
    const QR &qr;
  };

}

#endif //DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH
