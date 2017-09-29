// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>
#include <dune/common/hybridutilities.hh>

#include "../brezzidouglasmarini/brezzidouglasmarini1cube2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini1cube3d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini2cube2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh"
#include <dune/localfunctions/dualmortarbasis.hh>
#include "../refined/refinedp1.hh"
#include "../refined/refinedp0.hh"
#include "../hierarchical/hierarchicalp2.hh"
#include "../hierarchical/hierarchicalp2withelementbubble.hh"
#include "../hierarchical/hierarchicalprismp2.hh"
#include "../rannacherturek/rannacherturek.hh"
#include "../raviartthomas/raviartthomassimplex.hh"
#include "../raviartthomas/raviartthomascube.hh"

#include "test-localfe.hh"

int main(int argc, char** argv) try
{
  bool success = true;

  Dune::RefinedP1LocalFiniteElement<double,double,1> refp11dlfem;
  TEST_FE(refp11dlfem);

  Dune::RefinedP1LocalFiniteElement<double,double,2> refp12dlfem;
  TEST_FE(refp12dlfem);

  Dune::RefinedP1LocalFiniteElement<double,double,3> refp13dlfem;
  TEST_FE(refp13dlfem);

  Dune::RefinedP0LocalFiniteElement<double,double,1> refp01dlfem;
  TEST_FE(refp01dlfem);

  Dune::RefinedP0LocalFiniteElement<double,double,2> refp02dlfem;
  TEST_FE(refp02dlfem);

  // --------------------------------------------------------
  //  Test Brezzi-Douglas-Marini finite elements
  // --------------------------------------------------------
  Dune::BDM1Cube2DLocalFiniteElement<double,double> bdm1cube2dlfem(1);
  TEST_FE(bdm1cube2dlfem);

  Dune::BDM1Cube3DLocalFiniteElement<double,double> bdm1cube3dlfem(1);
  TEST_FE2(bdm1cube3dlfem, DisableLocalInterpolation);

  Dune::BDM2Cube2DLocalFiniteElement<double,double> bdm2cube2dlfem(1);
  TEST_FE(bdm2cube2dlfem);

  Dune::BDM1Simplex2DLocalFiniteElement<double,double> bdm1simplex2dlfem(1);
  TEST_FE(bdm1simplex2dlfem);

  Dune::BDM2Simplex2DLocalFiniteElement<double,double> bdm2simplex2dlfem(1);
  TEST_FE(bdm2simplex2dlfem);

  // --------------------------------------------------------
  //  Test hierarchical P2 finite elements
  // --------------------------------------------------------
  Dune::HierarchicalP2LocalFiniteElement<double,double,1> hierarchicalp21dlfem;
  TEST_FE(hierarchicalp21dlfem);

  Dune::HierarchicalP2LocalFiniteElement<double,double,2> hierarchicalp22dlfem;
  TEST_FE(hierarchicalp22dlfem);

  Dune::HierarchicalP2LocalFiniteElement<double,double,3> hierarchicalp23dlfem;
  TEST_FE(hierarchicalp23dlfem);

  Dune::HierarchicalPrismP2LocalFiniteElement<double,double> hierarchicalprismp2lfem;
  TEST_FE(hierarchicalprismp2lfem);

  Dune::HierarchicalP2WithElementBubbleLocalFiniteElement<double,double,2> hierarchicalp2bubble2dlfem;
  TEST_FE(hierarchicalp2bubble2dlfem);


  // --------------------------------------------------------
  //  Test dual mortar elements
  // --------------------------------------------------------
  Dune::DualP1LocalFiniteElement<double,double,1> dualp11dlfem;
  TEST_FE(dualp11dlfem);

  Dune::DualP1LocalFiniteElement<double,double,2> dualp12dlfem;
  TEST_FE(dualp12dlfem);

  Dune::DualP1LocalFiniteElement<double,double,3> dualp13dlfem;
  TEST_FE(dualp13dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,1> dualq11dlfem;
  TEST_FE(dualq11dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,2> dualq12dlfem;
  TEST_FE(dualq12dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,3> dualq13dlfem;
  TEST_FE(dualq13dlfem);


  // --------------------------------------------------------
  //  Test Raviart-Thomas Finite elements
  // --------------------------------------------------------
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

  // --------------------------------------------------------
  //  Test Rannacher-Turek Finite elements
  // --------------------------------------------------------
  Dune::RannacherTurekLocalFiniteElement<double,double,2> rannacher_turek2dfem;
  TEST_FE(rannacher_turek2dfem);

  Dune::RannacherTurekLocalFiniteElement<double,double,3> rannacher_turek3dfem;
  TEST_FE(rannacher_turek3dfem);

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
