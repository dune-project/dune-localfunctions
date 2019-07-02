// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <iostream>

#include <dune/localfunctions/raviartthomas/raviartthomassimplex.hh>
#include <dune/localfunctions/raviartthomas/raviartthomascube.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomaslfecache.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double> rt0simplex2dlfem(Dune::GeometryTypes::simplex(2),0);
  TEST_FE(rt0simplex2dlfem);

  Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double> rt1simplex2dlfem(Dune::GeometryTypes::simplex(2),1);
  TEST_FE(rt1simplex2dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,0> rt0cube2dlfem(1);
  TEST_FE(rt0cube2dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,0> rt0cube3dlfem(1);
  TEST_FE(rt0cube3dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,1> rt1cube2dlfem(1);
  TEST_FE(rt1cube2dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,1> rt1cube3dlfem(1);
  TEST_FE(rt1cube3dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,2> rt2cube2dlfem(1);
  TEST_FE(rt2cube2dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,3> rt3cube2dlfem(1);
  TEST_FE(rt3cube2dlfem);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,4> rt4cube2dlfem(1);
  TEST_FE(rt4cube2dlfem);

  Dune::RT0Cube2DLocalFiniteElement<double,double> rt0cube2dlfemDedicated;
  TEST_FE(rt0cube2dlfemDedicated);

  Dune::RT1Cube2DLocalFiniteElement<double,double> rt1cube2dlfemDedicated;
  TEST_FE(rt1cube2dlfemDedicated);

  Dune::RT02DLocalFiniteElement<double,double> rt02dlfemDedicated;
  TEST_FE(rt02dlfemDedicated);

  Dune::RT12DLocalFiniteElement<double,double> rt12dlfemDedicated;
  TEST_FE(rt12dlfemDedicated);

  // Test the RaviartThomasLocalFiniteElementCache
  Dune::RaviartThomasLocalFiniteElementCache<double,double,2,0> lagrangeLFECache;
  TEST_FE(lagrangeLFECache.get(Dune::GeometryTypes::simplex(2)));
  TEST_FE(lagrangeLFECache.get(Dune::GeometryTypes::cube(2)));

  // Test whether asking the cache for an element of the wrong dimension throws an exception
  bool lagrangeLFESuccess = false;
  try {
    auto doesntExist = lagrangeLFECache.get(Dune::GeometryTypes::simplex(1));
  } catch (Dune::Exception& e)
  {
    lagrangeLFESuccess = true;
  }
  success &= lagrangeLFESuccess;

  return success ? 0 : 1;
}
