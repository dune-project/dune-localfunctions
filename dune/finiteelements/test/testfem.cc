// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <typeinfo>
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
#include "../refinedp1.hh"

double epsilon = 1e-14;
double sqrt_epsilon = std::sqrt(epsilon);

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

template<class FE>
class LocalFEFunction
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType CT;

  LocalFEFunction(const FE fe) :
    fe_(fe)
  {
    resetCoefficients();
  }

  void resetCoefficients()
  {
    coeff_.resize(fe_.localBasis().size());
    for(int i=0; i<coeff_.size(); ++i)
      coeff_[i] = 0;
  }


  template<typename DT, typename RT>
  void evaluate (const DT& x, RT& y) const
  {
    typedef typename FE::Traits::LocalBasisType::Traits::RangeType RT;

    std::vector<RT> yy;
    fe_.localBasis().evaluateFunction(x, yy);

    y = 0.0;
    for (int i=0; i<yy.size(); ++i)
      y.axpy(coeff_[i], yy[i]);
  }

  std::vector<CT> coeff_;

private:
  const FE fe_;
};


template<class FE>
bool testLocalInterpolation(const FE& fe)
{
  bool success = true;
  LocalFEFunction<FE> f(fe);

  std::vector<typename LocalFEFunction<FE>::CT> coeff;
  for(int i=0; i<f.coeff_.size(); ++i)
  {
    f.coeff_[i] = 1;
    fe.localInterpolation().interpolate(f, coeff);
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type "
                << typeid(FE).name() << std::endl;
      std::cout << "    Interpolation vector has size " << coeff.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      success = false;
    }


    double diff=0;
    for(int j=0; j<coeff.size(); ++j)
    {
      if (std::abs(coeff[i]-f.coeff_[i]) > sqrt_epsilon)
      {
        std::cout << "Bug in LocalInterpolation for finite element type "
                  << typeid(FE).name() << std::endl;
        std::cout << "    Interpolation weight " << j
                  << " differs significantly from coefficient of linear combination." << std::endl;
        success = false;
      }
    }
  }
  return success;
}


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
  Dune::RefinedP1LocalFiniteElement<double,double> refp1lfem;

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
  refp1lfem.localInterpolation().interpolate(Func(),c);

  bool success = true;

  success = (testLocalInterpolation(p0lfem) and success);
  success = (testLocalInterpolation(p1lfem) and success);
  success = (testLocalInterpolation(p11dlfem) and success);
  success = (testLocalInterpolation(p12dlfem) and success);
  success = (testLocalInterpolation(p13dlfem) and success);
  success = (testLocalInterpolation(pk2dlfem) and success);
  success = (testLocalInterpolation(q1lfem) and success);
  success = (testLocalInterpolation(q12dlfem) and success);
  success = (testLocalInterpolation(q13dlfem) and success);
  success = (testLocalInterpolation(refp1lfem) and success);

  return success;
}
