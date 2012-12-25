// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_QKLOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  /** \todo Please doc me! */
  template<int k, int d, class LB>
  class QkLocalInterpolation
  {

    // Return i as a d-digit number in the (k+1)-nary system
    static Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(StaticPower<k+1,d>::power);

      for (int i=0; i<StaticPower<k+1,d>::power; i++)
      {
        // convert index i to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(i));

        // Generate coordinate of the i-th Lagrange point
        for (int j=0; j<d; j++)
          x[j] = (1.0*alpha[j])/k;

        f.evaluate(x,y); out[i] = y;
      }
    }
  };

  /** \todo Please doc me! */
  template<int d, class LB>
  class QkLocalInterpolation<0,d,LB>
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x(0);
      typename LB::Traits::RangeType y;
      f.evaluate(x,y);
      out.resize(1);
      out[0] = y;
    }
  };

}


#endif
