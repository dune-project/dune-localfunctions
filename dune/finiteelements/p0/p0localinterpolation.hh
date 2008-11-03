// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALINTERPOLATION_HH
#define DUNE_P0LOCALINTERPOLATION_HH

#include <dune/common/geometrytype.hh>

#include "../common/localinterpolation.hh"

namespace Dune
{

  template<class LB>
  class P0LocalInterpolation : LocalInterpolationInterface<P0LocalInterpolation>
  {
  public:

    //! determine coefficients interpolating a given function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typedef typename LB::Traits::DomainType DomainType;
      typedef typename LB::Traits::RangeType RangeType;
      typedef typename LB::Traits::DomainFieldType DF;
      const int dim=LB::Traits::dimDomain;

      DomainType x = Dune::ReferenceElements<DF,dim>::general(e.type()).position(0,0);
      RangeType y;

      out.resize(1);
      f.eval_local(x,y); out[0] = y;
    }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };

}

#endif
