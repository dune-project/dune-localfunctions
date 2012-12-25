// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASBASIS_HH
#define DUNE_RAVIARTTHOMASBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/utility/defaultbasisfactory.hh>
#include "raviartthomassimplexinterpolation.hh"
#include "raviartthomassimplexprebasis.hh"

namespace Dune
{
  template< unsigned int dim, class SF, class CF >
  struct RaviartThomasBasisFactory
    : public DefaultBasisFactory< RTPreBasisFactory<dim,CF>,
          RaviartThomasL2InterpolationFactory<dim,CF>,
          dim,dim,SF,CF >
  {};
}

#endif // #ifndef DUNE_RAVIARTTHOMASBASIS_HH
