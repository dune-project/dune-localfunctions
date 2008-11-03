// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>

#include "../p0.hh"
#include "../p12d.hh"
#include "../pk2d.hh"
#include "../q12d.hh"
#include "../q22d.hh"

class Func
{
public:
  template<typename DT, typename RT>
  void evaluate (const DT& x, RT& y) const
  {
    DT c(0.5);
    c -= x;
    y[0] = exp(-3.0*c.two_norm2());
  }
};

int main(int argc, char** argv)
{
  Dune::P0LocalFiniteElement<double,double,2> p0lfem(Dune::GeometryType::simplex);
  Dune::P12DLocalFiniteElement<double,double> p12dlfem;
  Dune::Pk2DLocalFiniteElement<double,double,5> pk2dlfem(3);
  Dune::Q12DLocalFiniteElement<double,double> q12dlfem;
  Dune::Q22DLocalFiniteElement<double,double> q22dlfem;

  std::vector<double> c;

  p0lfem.localInterpolation().interpolate(Func(),c);
  p12dlfem.localInterpolation().interpolate(Func(),c);
  pk2dlfem.localInterpolation().interpolate(Func(),c);
  q12dlfem.localInterpolation().interpolate(Func(),c);
  q22dlfem.localInterpolation().interpolate(Func(),c);

  return 0;
}
