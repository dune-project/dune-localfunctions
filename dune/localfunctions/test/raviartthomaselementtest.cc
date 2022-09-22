// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include "config.h"

#include <iostream>

#include <dune/localfunctions/raviartthomas/raviartthomassimplex.hh>
#include <dune/localfunctions/raviartthomas/raviartthomascube.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas03d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomaslfecache.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double> rt0simplex2dlfem(Dune::GeometryTypes::simplex(2),0);
  TEST_FE3(rt0simplex2dlfem, DisableNone, 2);

  Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double> rt1simplex2dlfem(Dune::GeometryTypes::simplex(2),1);
  TEST_FE3(rt1simplex2dlfem, DisableNone,2);

  Dune::RaviartThomasSimplexLocalFiniteElement<3,double,double> rt0simplex3dlfem(Dune::GeometryTypes::simplex(3),0);
  TEST_FE3(rt0simplex3dlfem,DisableNone,2);

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,0> rt0cube2dlfem;
  TEST_FE(rt0cube2dlfem);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,0> rt0cube2dlfem(s);
    TEST_FE(rt0cube2dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,0> rt0cube3dlfem;
  TEST_FE(rt0cube3dlfem);
  for (unsigned int s = 0; s < 64; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,0> rt0cube3dlfem(s);
    TEST_FE(rt0cube3dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,1> rt1cube2dlfem;
  TEST_FE(rt1cube2dlfem);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,1> rt1cube2dlfem(s);
    TEST_FE(rt1cube2dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,1> rt1cube3dlfem;
  TEST_FE(rt1cube3dlfem);
  for (unsigned int s = 0; s < 64; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,1> rt1cube3dlfem(s);
    TEST_FE(rt1cube3dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,2> rt2cube2dlfem;
  TEST_FE(rt2cube2dlfem);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,2> rt2cube2dlfem(s);
    TEST_FE(rt2cube2dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,3> rt3cube2dlfem;
  TEST_FE(rt3cube2dlfem);
  for (unsigned int s = 0; s < 64; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,3> rt3cube2dlfem(s);
    TEST_FE(rt3cube2dlfem);
  }

  Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,4> rt4cube2dlfem;
  TEST_FE(rt4cube2dlfem);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,4> rt4cube2dlfem(s);
    TEST_FE(rt4cube2dlfem);
  }

  Dune::RT0Cube2DLocalFiniteElement<double,double> rt0cube2dlfemDedicated;
  TEST_FE(rt0cube2dlfemDedicated);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RT0Cube2DLocalFiniteElement<double,double> rt0cube2dlfemDedicated(s);
    TEST_FE(rt0cube2dlfemDedicated);
  }

  Dune::RT1Cube2DLocalFiniteElement<double,double> rt1cube2dlfemDedicated;
  TEST_FE(rt1cube2dlfemDedicated);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RT1Cube2DLocalFiniteElement<double,double> rt1cube2dlfemDedicated(s);
    TEST_FE(rt1cube2dlfemDedicated);
  }

  Dune::RT02DLocalFiniteElement<double,double> rt02dlfemDedicated;
  TEST_FE(rt02dlfemDedicated);
  for (unsigned int s = 0; s < 8; s++)
  {
    Dune::RT02DLocalFiniteElement<double,double> rt02dlfemDedicated(s);
    TEST_FE(rt02dlfemDedicated);
  }

  Dune::RT12DLocalFiniteElement<double,double> rt12dlfemDedicated;
  TEST_FE(rt12dlfemDedicated);
  for (unsigned int s = 0; s < 8; s++)
  {
    Dune::RT12DLocalFiniteElement<double,double> rt12dlfemDedicated(s);
    TEST_FE(rt12dlfemDedicated);
  }

  Dune::RT03DLocalFiniteElement<double,double> rt03dlfemDedicated;
  TEST_FE(rt03dlfemDedicated);
  for (unsigned int s = 0; s < 16; s++)
  {
    Dune::RT03DLocalFiniteElement<double,double> rt03dlfemDedicated(s);
    TEST_FE(rt03dlfemDedicated);
  }

  Dune::RT0Cube3DLocalFiniteElement<double,double> rt0cube3dlfemDedicated;
  TEST_FE(rt0cube3dlfemDedicated);
  for (unsigned int s = 0; s < 64; s++)
  {
    Dune::RT0Cube3DLocalFiniteElement<double,double> rt0cube3dlfemDedicated(s);
    TEST_FE(rt0cube3dlfemDedicated);
  }

  Dune::RT1Cube3DLocalFiniteElement<double,double> rt1cube3dlfemDedicated;
  TEST_FE(rt1cube3dlfemDedicated);
  for (unsigned int s = 0; s < 64; s++)
  {
    Dune::RT1Cube3DLocalFiniteElement<double,double> rt1cube3dlfemDedicated(s);
    TEST_FE(rt1cube3dlfemDedicated);
  }

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
