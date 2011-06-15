// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q2_LOCALINTERPOLATION_HH
#define DUNE_Q2_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class Q2LocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      static const int dim = LB::Traits::dimDomain;

      switch (dim) {

      case 3 : {

        out.resize(27);
        x[0] = 0.0; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y); out[0] = y;
        x[0] = 0.5; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y); out[14] = y;
        x[0] = 1.0; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y); out[1] = y;
        x[0] = 0.0; x[1] = 0.5; x[2] = 0.0; f.evaluate(x,y); out[12] = y;
        x[0] = 0.5; x[1] = 0.5; x[2] = 0.0; f.evaluate(x,y); out[24] = y;
        x[0] = 1.0; x[1] = 0.5; x[2] = 0.0; f.evaluate(x,y); out[13] = y;
        x[0] = 0.0; x[1] = 1.0; x[2] = 0.0; f.evaluate(x,y); out[2] = y;
        x[0] = 0.5; x[1] = 1.0; x[2] = 0.0; f.evaluate(x,y); out[15] = y;
        x[0] = 1.0; x[1] = 1.0; x[2] = 0.0; f.evaluate(x,y); out[3] = y;

        x[0] = 0.0; x[1] = 0.0; x[2] = 0.5; f.evaluate(x,y); out[8] = y;
        x[0] = 0.5; x[1] = 0.0; x[2] = 0.5; f.evaluate(x,y); out[22] = y;
        x[0] = 1.0; x[1] = 0.0; x[2] = 0.5; f.evaluate(x,y); out[9] = y;
        x[0] = 0.0; x[1] = 0.5; x[2] = 0.5; f.evaluate(x,y); out[20] = y;
        x[0] = 0.5; x[1] = 0.5; x[2] = 0.5; f.evaluate(x,y); out[26] = y;
        x[0] = 1.0; x[1] = 0.5; x[2] = 0.5; f.evaluate(x,y); out[21] = y;
        x[0] = 0.0; x[1] = 1.0; x[2] = 0.5; f.evaluate(x,y); out[10] = y;
        x[0] = 0.5; x[1] = 1.0; x[2] = 0.5; f.evaluate(x,y); out[23] = y;
        x[0] = 1.0; x[1] = 1.0; x[2] = 0.5; f.evaluate(x,y); out[11] = y;

        x[0] = 0.0; x[1] = 0.0; x[2] = 1.0; f.evaluate(x,y); out[4] = y;
        x[0] = 0.5; x[1] = 0.0; x[2] = 1.0; f.evaluate(x,y); out[18] = y;
        x[0] = 1.0; x[1] = 0.0; x[2] = 1.0; f.evaluate(x,y); out[5] = y;
        x[0] = 0.0; x[1] = 0.5; x[2] = 1.0; f.evaluate(x,y); out[16] = y;
        x[0] = 0.5; x[1] = 0.5; x[2] = 1.0; f.evaluate(x,y); out[25] = y;
        x[0] = 1.0; x[1] = 0.5; x[2] = 1.0; f.evaluate(x,y); out[17] = y;
        x[0] = 0.0; x[1] = 1.0; x[2] = 1.0; f.evaluate(x,y); out[6] = y;
        x[0] = 0.5; x[1] = 1.0; x[2] = 1.0; f.evaluate(x,y); out[19] = y;
        x[0] = 1.0; x[1] = 1.0; x[2] = 1.0; f.evaluate(x,y); out[7] = y;

        break;
      }
      default :
        DUNE_THROW(NotImplemented, "Q2LocalInterpolation for dim==" << dim);
      }
    }
  };
}

#endif
