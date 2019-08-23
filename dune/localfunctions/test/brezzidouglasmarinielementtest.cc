// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <iostream>

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarinicube.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,2,1> bdm1cube2dlfem(1);
  TEST_FE(bdm1cube2dlfem);

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,3,1> bdm1cube3dlfem(1);
  TEST_FE2(bdm1cube3dlfem, DisableLocalInterpolation);

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,2,2> bdm2cube2dlfem(1);
  TEST_FE(bdm2cube2dlfem);

  Dune::BDM1Simplex2DLocalFiniteElement<double,double> bdm1simplex2dlfem(1);
  TEST_FE(bdm1simplex2dlfem);

  Dune::BDM2Simplex2DLocalFiniteElement<double,double> bdm2simplex2dlfem(1);
  TEST_FE(bdm2simplex2dlfem);

  return success ? 0 : 1;
}
