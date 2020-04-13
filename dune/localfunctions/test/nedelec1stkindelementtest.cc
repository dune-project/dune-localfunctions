// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/localfunctions/nedelec/nedelec1stkindsimplex.hh>

#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;

int main(int argc, char** argv)
{
  bool success = true;

  {
    Nedelec1stKindSimplexLocalFiniteElement<double,double,2,1> nedelecLFEM;
    TEST_FE(nedelecLFEM);
  }

  return success ? 0 : 1;
}
