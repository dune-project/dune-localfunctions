// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_PRISM_P2_LOCALINTERPOLATION_HH
#define DUNE_HIERARCHICAL_PRISM_P2_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  /**
     \tparam LB The LocalBasis implementation
   */
  template<class LB>
  class HierarchicalPrismP2LocalInterpolation
  {
  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      out.resize(18);

      //First the  vertex dofs
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.0;    f.evaluate(x, y);    out[0] = y;
      x[0] = 1.0;    x[1] = 0.0;     x[2] = 0.0;    f.evaluate(x, y);    out[1] = y;
      x[0] = 0.0;    x[1] = 1.0;     x[2] = 0.0;    f.evaluate(x, y);    out[2] = y;
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 1.0;    f.evaluate(x, y);    out[3] = y;
      x[0] = 1.0;    x[1] = 0.0;     x[2] = 1.0;    f.evaluate(x, y);    out[4] = y;
      x[0] = 0.0;    x[1] = 1.0;     x[2] = 1.0;    f.evaluate(x, y);    out[5] = y;


      // Then: the 9 edge dofs and the 3 face dofs
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.5;    f.evaluate(x, y);
      out[6] = y - 0.5*(out[0] + out[3]);

      x[0] = 1.0;    x[1] = 0.0;     x[2] = 0.5;    f.evaluate(x, y);
      out[7] = y - 0.5*(out[1] + out[4]);

      x[0] = 0.0;    x[1] = 1.0;     x[2] = 0.5;    f.evaluate(x, y);
      out[8] = y - 0.5*(out[2] + out[5]);

      x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.0;    f.evaluate(x, y);
      out[9] = y - 0.5*(out[0] + out[1]);

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.0;    f.evaluate(x, y);
      out[10] = y - 0.5*(out[2] + out[0]);

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 0.0;    f.evaluate(x, y);
      out[11] = y - 0.5*(out[2] + out[1]);

      x[0] = 0.5;    x[1] = 0.0;     x[2] = 1.0;    f.evaluate(x, y);
      out[12] = y - 0.5*(out[3] + out[4]);

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 1.0;    f.evaluate(x, y);
      out[13] = y - 0.5*(out[3] + out[5]);

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 1.0;    f.evaluate(x, y);
      out[14] = y - 0.5*(out[4] + out[5]);


      //faces
      x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.5;    f.evaluate(x, y);
      out[15] = y - 0.25*(out[4] + out[1] + out[0] +  out[3] );

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.5;    f.evaluate(x, y);
      out[16] = y - 0.25*(out[2] + out[0] + out[3] + out[5] );

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 0.5;    f.evaluate(x, y);
      out[17] = y - 0.25*(out[2] + out[1] + out[4] + out[5] );

    }
  };
}

#endif
