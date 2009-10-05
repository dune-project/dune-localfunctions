// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "orthonormalbasis.hh"
using namespace Dune;
using namespace GenericGeometry;
int main ( int argc, char **argv )
{
  int order = atoi(argv[1]);
  OrthonormalBasis<TOPOLOGY,double> onb(order);
}
