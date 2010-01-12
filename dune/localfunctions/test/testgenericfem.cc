// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>

#include <dune/localfunctions/generic/math/gmpfield.hh>

#include <dune/localfunctions/generic/lagrangefiniteelement.hh>
#include <dune/localfunctions/generic/lagrangebasis/equidistantpoints.hh>

#include "testfem.hh"


int main(int argc, char** argv) try
{
  bool success = true;

  for (unsigned int order=1; order<=4; ++order)
  {
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double>
    lagrangeSimplex(0b00,order);
    success &= testFE(lagrangeSimplex);
  }
  for (unsigned int order=1; order<=3; ++order)
  {
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double>
    lagrangeCube(0b11,order);
    success &= testFE(lagrangeCube);
  }
#if HAVE_GMP
  for (unsigned int order=5; order<=10; ++order)
  {
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double,
        Dune::GMPField<128>,Dune::GMPField<256> >
    lagrangeSimplex(0b00,order);
    success &= testFE(lagrangeSimplex);
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double,
        Dune::GMPField<128>,Dune::GMPField<256> >
    lagrangeCube(0b11,order);
    success &= testFE(lagrangeCube);
  }
#endif

  return success;
}
catch (const Dune::Exception &e)
{
  std::cout << e << std::endl;
  return 1;
}
