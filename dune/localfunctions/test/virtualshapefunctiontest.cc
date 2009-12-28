// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstddef>
#include <iostream>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/localfunctions/p1.hh>

/** \file
    \brief Test the dynamically polymorphic shape function interface
 */

using namespace Dune;


int main (int argc, char *argv[]) try
{

  const Dune::LocalFiniteElementVirtualImp<Dune::P1LocalFiniteElement<double, double, 2> > virtualLocalSimplexFE_;

  return 0;

}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
