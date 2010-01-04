// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//#define DUNE_VIRTUAL_SHAPEFUNCTIONS 1

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>

#include "../p0.hh"
#include "../p1.hh"
#include "../prismp1.hh"
#include "../prismp2.hh"
#include "../q1.hh"
#include "../refinedp1.hh"
#include "../refinedp0.hh"
#include "../p23d.hh"
#include "../hierarchicalp2.hh"
#include "../hierarchicalp2withelementbubble.hh"
#include "../hierarchicalprismp2.hh"
#include "../rannacher_turek2d.hh"

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
// these shape functions don't provide the virtual interface
#include "../pk2d.hh"
#include "../pk3d.hh"
#include "../q22d.hh"
#include "../rt02d.hh"
#include "../rt0q2d.hh"
#include "../rt0q3d.hh"
#include "../monom.hh"
#include "../edger02d.hh"
#include "../edges02d.hh"
#include "../edges03d.hh"
#endif


double TOL = 1e-10;

template<class FE>
class Func :
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  public Dune::Function<
      const typename FE::Traits::LocalBasisType::Traits::DomainType&,
      typename FE::Traits::LocalBasisType::Traits::RangeType&>
#else
  public Dune::VirtualFunction<
      typename FE::Traits::LocalBasisType::Traits::DomainType,
      typename FE::Traits::LocalBasisType::Traits::RangeType>
#endif
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;

  // the base class exports the type but some FEs assume that
  // it is encapsulated in a Traits class
  struct Traits
  {
    typedef typename Base::RangeType RangeType;
    typedef typename Base::DomainType DomainType;
  };

  void evaluate (const DomainType& x, RangeType& y) const
  {
    DomainType c(0.5);

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
class LocalFEFunction :
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  public Dune::Function<
      const typename FE::Traits::LocalBasisType::Traits::DomainType&,
      typename FE::Traits::LocalBasisType::Traits::RangeType&>
#else
  public Dune::VirtualFunction<
      typename FE::Traits::LocalBasisType::Traits::DomainType,
      typename FE::Traits::LocalBasisType::Traits::RangeType>
#endif
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;

  // the base class exports the type but some FEs assume that
  // it is encapsulated in a Traits class
  struct Traits
  {
    typedef typename Base::RangeType RangeType;
    typedef typename Base::DomainType DomainType;
  };

  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType CT;

  LocalFEFunction(const FE& fe) :
    fe_(fe)
  {
    resetCoefficients();
  }

  void resetCoefficients()
  {
    coeff_.resize(fe_.localBasis().size());
    for(std::size_t i=0; i<coeff_.size(); ++i)
      coeff_[i] = 0;
  }

  void setRandom(double max)
  {
    coeff_.resize(fe_.localBasis().size());
    for(std::size_t i=0; i<coeff_.size(); ++i)
      coeff_[i] = ((1.0*std::rand()) / RAND_MAX - 0.5)*2.0*max;
  }


  void evaluate (const DomainType& x, RangeType& y) const
  {
    std::vector<RangeType> yy;
    fe_.localBasis().evaluateFunction(x, yy);

    y = 0.0;
    for (std::size_t i=0; i<yy.size(); ++i)
      y.axpy(coeff_[i], yy[i]);
  }

  std::vector<CT> coeff_;

private:
  const FE& fe_;
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
    for(std::size_t j=0; j<coeff.size(); ++j)
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

  fe.localInterpolation().interpolate(Func<FE>(),c);

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  return testLocalInterpolation<FE>(fe);
#else
  typedef typename FE::Traits::LocalBasisType::Traits::DomainFieldType DT;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType RT;
  const int dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  typedef Dune::LocalFiniteElementInterface<DT, RT, dim> FEBase;
  return testLocalInterpolation<FEBase>(fe);
#endif
}


// tmp for testing arbitrary order finite elements
template<int k>
bool testArbitraryOrderFE()
{
  bool success = true;

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  Dune::Pk2DLocalFiniteElement<double,double,k> pk2dlfem(1);
  success = testFE(pk2dlfem) and success;

  Dune::Pk3DLocalFiniteElement<double,double,k> pk3dlfem;
  success = testFE(pk3dlfem) and success;
#endif

  return testArbitraryOrderFE<k-1>() and success;
}

template<>
bool testArbitraryOrderFE<0>()
{
  return true;
}

template<int k>
bool testMonomials()
{
  bool success = true;

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  Dune::MonomLocalFiniteElement<double,double,1,k> monom1d(Dune::GeometryType::simplex);
  success = testFE(monom1d) and success;

  Dune::MonomLocalFiniteElement<double,double,2,k> monom2d(Dune::GeometryType::simplex);
  success = testFE(monom2d) and success;

  Dune::MonomLocalFiniteElement<double,double,3,k> monom3d(Dune::GeometryType::simplex);
  success = testFE(monom3d) and success;
#endif

  return testMonomials<k-1>() and success;
}

template<>
bool testMonomials<-1>()
{
  return true;
}

int main(int argc, char** argv) try
{
  bool success = true;

  Dune::P0LocalFiniteElement<double,double,2> p0lfem(Dune::GeometryType::simplex);
  success = testFE(p0lfem) and success;

  Dune::P1LocalFiniteElement<double,double,1> p11dlfem;
  success = testFE(p11dlfem) and success;

  Dune::P1LocalFiniteElement<double,double,2> p12dlfem;
  success = testFE(p12dlfem) and success;

  Dune::P1LocalFiniteElement<double,double,3> p13dlfem;
  success = testFE(p13dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,1> q11dlfem;
  success = testFE(q11dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,2> q12dlfem;
  success = testFE(q12dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,3> q13dlfem;
  success = testFE(q13dlfem) and success;

  Dune::RefinedP1LocalFiniteElement<double,double,2> refp12dlfem;
  success = testFE(refp12dlfem) and success;

  Dune::RefinedP1LocalFiniteElement<double,double,3> refp13dlfem;
  success = testFE(refp13dlfem) and success;

  Dune::RefinedP0LocalFiniteElement<double,double,2> refp02dlfem;
  success = testFE(refp02dlfem) and success;

  Dune::P23DLocalFiniteElement<double,double> p23dlfem;
  success = testFE(p23dlfem) and success;

  //    Dune::HierarchicalP2LocalFiniteElement<double,double,1> hierarchicalp21dlfem;
  //    success = testFE(hierarchicalp21dlfem) and success;

  Dune::HierarchicalP2LocalFiniteElement<double,double,2> hierarchicalp22dlfem;
  success = testFE(hierarchicalp22dlfem) and success;

  Dune::HierarchicalP2LocalFiniteElement<double,double,3> hierarchicalp23dlfem;
  success = testFE(hierarchicalp23dlfem) and success;

  Dune::HierarchicalPrismP2LocalFiniteElement<double,double> hierarchicalprismp2lfem;
  success = testFE(hierarchicalprismp2lfem) and success;

  Dune::HierarchicalP2WithElementBubbleLocalFiniteElement<double,double,2> hierarchicalp2bubble2dlfem;
  success = testFE(hierarchicalp2bubble2dlfem) and success;

  Dune::PrismP1LocalFiniteElement<double,double> prismp1fem;
  success = testFE(prismp1fem) and success;

  Dune::PrismP2LocalFiniteElement<double,double> prismp2fem;
  success = testFE(prismp2fem) and success;

  success = testArbitraryOrderFE<10>() and success;

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  Dune::Q22DLocalFiniteElement<double,double> q22dlfem;
  success = testFE(q22dlfem) and success;

  Dune::EdgeR02DLocalFiniteElement<double,double> edger02dlfem;
  success = testFE(edger02dlfem) and success;

  Dune::EdgeS02DLocalFiniteElement<double,double> edges02dlfem;
  success = testFE(edges02dlfem) and success;

  Dune::EdgeS03DLocalFiniteElement<double,double> edges03dlfem;
  success = testFE(edges03dlfem) and success;

  Dune::RT02DLocalFiniteElement<double,double> rt02dlfem(1);
  success = testFE(rt02dlfem) and success;

  Dune::RT0Q2DLocalFiniteElement<double,double> rt0q2dlfem(1);
  success = testFE(rt0q2dlfem) and success;

  Dune::RT0Q3DLocalFiniteElement<double,double> rt0q3dlfem(1);
  success = testFE(rt0q3dlfem) and success;

  Dune::RannacherTurek2DLocalFiniteElement<double,double> rannacher_turek2dfem;
  success = testFE(rannacher_turek2dfem) and success;

#endif

  std::cout << "Monomials are only tested up to order 2 due to the instability of interpolate()." << std::endl;
  success = testMonomials<2>() and success;


  return success;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
