// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>

#include "../p0.hh"
#include "../p1.hh"
#include "../p11d.hh"
#include "../p12d.hh"
#include "../p13d.hh"
#include "../pk2d.hh"
#include "../q1.hh"
#include "../q12d.hh"
#include "../q13d.hh"
#include "../q22d.hh"
#include "../rt02d.hh"

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
  Dune::P1LocalFiniteElement<double,double,2> p1lfem;
  Dune::P12DLocalFiniteElement<double,double> p11dlfem;
  Dune::P12DLocalFiniteElement<double,double> p12dlfem;
  Dune::P13DLocalFiniteElement<double,double> p13dlfem;
  Dune::Pk2DLocalFiniteElement<double,double,5> pk2dlfem(3);
  Dune::Q1LocalFiniteElement<double,double,3> q1lfem;
  Dune::Q12DLocalFiniteElement<double,double> q12dlfem;
  Dune::Q12DLocalFiniteElement<double,double> q13dlfem;
  Dune::Q22DLocalFiniteElement<double,double> q22dlfem;
  Dune::RT02DLocalFiniteElement<double,double> rt02dlfem;

  std::vector<double> c;

  p0lfem.localInterpolation().interpolate(Func(),c);
  p1lfem.localInterpolation().interpolate(Func(),c);
  p11dlfem.localInterpolation().interpolate(Func(),c);
  p12dlfem.localInterpolation().interpolate(Func(),c);
  p13dlfem.localInterpolation().interpolate(Func(),c);
  pk2dlfem.localInterpolation().interpolate(Func(),c);
  q1lfem.localInterpolation().interpolate(Func(),c);
  q12dlfem.localInterpolation().interpolate(Func(),c);
  q13dlfem.localInterpolation().interpolate(Func(),c);
  q22dlfem.localInterpolation().interpolate(Func(),c);

  return 0;
}
