// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/utility/defaultbasisfactory.hh>
#include <dune/localfunctions/utility/monomialbasis.hh>

#include <dune/localfunctions/lagrange/interpolation.hh>

namespace Dune
{

  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisFactory
    : public DefaultBasisFactory< MonomialBasisFactory<dim,CF>,
          LagrangeInterpolationFactory<LP,dim,CF>,
          dim,1,SF,CF >
  {};

}

#endif // #ifndef DUNE_LAGRANGEBASIS_HH
