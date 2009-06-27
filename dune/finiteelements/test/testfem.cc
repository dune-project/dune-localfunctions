// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include "../p0.hh"
#include "../p1.hh"
#include "../p11d.hh"
#include "../p12d.hh"
#include "../p13d.hh"
#include "../pk2d.hh"
#include "../p23d.hh"
#include "../pk3d.hh"
#include "../q1.hh"
#include "../q12d.hh"
#include "../q13d.hh"
#include "../q22d.hh"
#include "../rt02d.hh"
#include "../refinedp1.hh"
#include "../monom.hh"


double TOL = 1e-4;

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

// This class defines a local finite element function.
// It is determined by a local finite element and
// representing the local basis and a coefficient vector.
// This provides the evaluate method needed by the interpolate()
// method.
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

  void setRandom(double max)
  {
    coeff_.resize(fe_.localBasis().size());
    for(int i=0; i<coeff_.size(); ++i)
      coeff_[i] = ((std::rand() / RAND_MAX) - 0.5)*2.0*max;
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


// Check if localInterpolation is consistens with
// localBasis evaluation.
template<class FE>
bool testLocalInterpolation(const FE& fe, int n=100)
{
  bool success = true;
  LocalFEFunction<FE> f(fe);

  std::vector<typename LocalFEFunction<FE>::CT> coeff;
  for(int i=0; i<n; ++i)
  {
    // Set random coefficient vector
    f.setRandom(100);

    // Compute interpolation weights
    fe.localInterpolation().interpolate(f, coeff);

    // Check size of weight vector
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type "
                << typeid(FE).name() << std::endl;
      std::cout << "    Interpolation vector has size " << coeff.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      success = false;
    }

    // Check if interpolation weights are equal to coefficients
    for(int j=0; j<coeff.size(); ++j)
    {
      if (std::abs(coeff[j]-f.coeff_[j]) > TOL)
      {
        std::cout << "Bug in LocalInterpolation for finite element type "
                  << typeid(FE).name() << std::endl;
        std::cout << "    Interpolation weight " << j
                  << " differs by " << std::abs(coeff[j]-f.coeff_[j])
                  << " from coefficient of linear combination." << std::endl;
        success = false;
      }
    }
  }
  return success;
}


// call tests for given finite element
template<class FE>
bool testFE(const FE& fe)
{
  std::vector<double> c;
  fe.localInterpolation().interpolate(Func(),c);

  return testLocalInterpolation(fe);
}


// tmp for testing arbitrary order finite elements
template<int k>
bool testArbitraryOrderFE()
{
  bool success = true;
  std::vector<double> c;

  Dune::Pk2DLocalFiniteElement<double,double,k> pk2dlfem(1);
  success = testFE(pk2dlfem) and success;

  Dune::Pk3DLocalFiniteElement<double,double,1> pk3dlfem;
  success = testFE(pk3dlfem) and success;

  Dune::MonomLocalFiniteElement<double,double,1,k> monom1d(Dune::GeometryType::simplex);
  success = testFE(monom1d) and success;

  Dune::MonomLocalFiniteElement<double,double,2,k> monom2d(Dune::GeometryType::simplex);
  success = testFE(monom2d) and success;

  Dune::MonomLocalFiniteElement<double,double,3,k> monom3d(Dune::GeometryType::simplex);
  success = testFE(monom3d) and success;

  return testArbitraryOrderFE<k-1>() and success;
}

template<>
bool testArbitraryOrderFE<0>()
{
  return true;
}

int main(int argc, char** argv)
{
  bool success = true;

  Dune::P0LocalFiniteElement<double,double,2> p0lfem(Dune::GeometryType::simplex);
  success = testFE(p0lfem) and success;

  Dune::P1LocalFiniteElement<double,double,2> p1lfem;
  success = testFE(p1lfem) and success;

  Dune::P12DLocalFiniteElement<double,double> p11dlfem;
  success = testFE(p11dlfem) and success;

  Dune::P12DLocalFiniteElement<double,double> p12dlfem;
  success = testFE(p12dlfem) and success;

  Dune::P13DLocalFiniteElement<double,double> p13dlfem;
  success = testFE(p13dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,3> q1lfem;
  success = testFE(q1lfem) and success;

  Dune::Q12DLocalFiniteElement<double,double> q12dlfem;
  success = testFE(q12dlfem) and success;

  Dune::Q12DLocalFiniteElement<double,double> q13dlfem;
  success = testFE(q13dlfem) and success;

  Dune::Q22DLocalFiniteElement<double,double> q22dlfem;
  success = testFE(q22dlfem) and success;

  Dune::RefinedP1LocalFiniteElement<double,double> refp1lfem;
  success = testFE(refp1lfem) and success;

  Dune::P23DLocalFiniteElement<double,double> p23dlfem;
  success = testFE(p23dlfem) and success;

  // Monomials produce an error for higher order since
  // the tolerance of 1e-4 is to small.
  success = testArbitraryOrderFE<5>() and success;

  Dune::RT02DLocalFiniteElement<double,double> rt02dlfem;

  return success;
}
