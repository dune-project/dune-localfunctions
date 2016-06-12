// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <typeinfo>
#include <fenv.h>

#include <dune/localfunctions/lagrange/p1.hh>
#include <dune/localfunctions/lagrange/p23d.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>
#include <dune/localfunctions/lagrange/pk3d.hh>

/** \file
    \brief Performs some tests for the Pk shape functions
 */

bool success = true;
double epsilon = 1e-14;
double sqrt_epsilon = std::sqrt(epsilon);

using namespace Dune;

// Generate list of Lagrange points for dim-dimensional simplex
template <int dim>
void getPkTestPoints(unsigned order, unsigned level, std::vector<FieldVector<double,dim> >& test_points)
{
  for (unsigned i = 0; i <= order - level; ++i)
  {
    std::vector<FieldVector<double,dim-1> > test_points_lower_dim;
    getPkTestPoints(order, level + i, test_points_lower_dim);
    double coord = double(i) / order;
    for (unsigned j = 0; j < test_points_lower_dim.size(); ++j)
    {
      FieldVector<double,dim> pos;
      for (int k = 0; k < dim-1; ++k)
        pos[k] = test_points_lower_dim[j][k];
      pos[dim-1] = coord;
      test_points.push_back(pos);
    }
  }
}

// Template specialization to terminate recursion
template <>
void getPkTestPoints(unsigned order, unsigned level, std::vector<FieldVector<double,0> >& test_points)
{
  FieldVector<double,0> pos;
  test_points.push_back(pos);
}

template <class FE>
void testPk(const FE& local_fe)
{
  typedef typename FE::Traits::LocalBasisType::Traits LBTraits;
  const int dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  const unsigned order = local_fe.localBasis().order();

  std::vector<FieldVector<double,1> > values;

  std::vector<FieldVector<double,dim> > test_points;
  getPkTestPoints(order, 0, test_points);

  for (unsigned n = 0; n < test_points.size(); ++n)
  {
    FieldVector<double,dim> pos = test_points[n];

    //////////////////////////////////////////////////////////////////
    //  Verfiy that shape functions fulfill \phi_i(x_j) = \delta_{ij}
    //  We assume that the shape functions are ordered corresponding
    //  to the test points returned by getPkTestPoints()
    //////////////////////////////////////////////////////////////////

    local_fe.localBasis().evaluateFunction(pos, values);
    for (unsigned i = 0; i < values.size(); ++i)
      if (std::abs(values[i] - double(i==n)) > epsilon)
      {
        std::cerr << "Bug in shape function in local finite element type"
                  << typeid(FE).name() << std::endl;
        std::cerr << "Shape function " << n << " has value " << values[i]
                  << " at position " << pos << " while " << double(i==n)
                  << " was expected" << std::endl;
        success = false;
      }

    //////////////////////////////////////////////////////////////////
    //  Check the partial derivatives by comparing them
    //  to finite difference approximations
    //////////////////////////////////////////////////////////////////

    // Get the shape function derivatives at pos
    std::vector<typename LBTraits::JacobianType> jacobians;
    local_fe.localBasis().evaluateJacobian(pos, jacobians);

    // Loop over all axes
    for (int k=0; k<dim; k++)
    {
      // Compute an approximation to the derivative by finite differences
      FieldVector<double,dim> upPos   = pos;
      FieldVector<double,dim> downPos = pos;

      upPos[k]   += sqrt_epsilon;
      downPos[k] -= sqrt_epsilon;

      std::vector<FieldVector<double,1> > upValues, downValues;

      local_fe.localBasis().evaluateFunction(upPos,   upValues);
      local_fe.localBasis().evaluateFunction(downPos, downValues);

      // Loop over all shape functions in this set
      for (unsigned j=0; j<local_fe.localBasis().size(); ++j)
      {
        // The current partial derivative, just for ease of notation
        double derivative = jacobians[j][0][k];

        // Compute finite difference approximation
        double finiteDiff = (upValues[j] - downValues[j]) / (2*sqrt_epsilon);

        // Check
        if (std::abs(derivative - finiteDiff) > sqrt_epsilon) {
          std::cerr << "Bug in shape function in local finite element type "
                    << typeid(FE).name() << std::endl;
          std::cerr << "    of order " << order << "." << std::endl;
          std::cerr << "    Shape function derivative differs "
                    << "significantly from FD approximation" << std::endl;
          std::cerr << "    Shape function " << j << " at position " << pos
                    << ":  derivative in direction " << k
                    << " is " << derivative << ", but " << finiteDiff
                    << " is expected." << std::endl;
          success = false;
        }
      }
    }
  }
}

int main (int argc, char *argv[]) try
{
#if __linux__
#if (!defined __INTEL_COMPILER || __INTEL_COMPILER >= 1010)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
#endif

  P1LocalFiniteElement<double,double,1> p11d;
  testPk(p11d);

  P1LocalFiniteElement<double,double,2> p12d;
  testPk(p12d);

  P1LocalFiniteElement<double,double,3> p13d;
  testPk(p13d);

  //     P23DLocalFiniteElement does not fulfill above assumption on the
  //     ordering of the shape functions
  //     P23DLocalFiniteElement<double,double> p23d;
  //     testPk(p23d);

  Pk2DLocalFiniteElement<double,double,1> pk12d;
  testPk(pk12d);
  Pk2DLocalFiniteElement<double,double,2> pk22d;
  testPk(pk22d);
  Pk2DLocalFiniteElement<double,double,3> pk32d;
  testPk(pk32d);
  Pk2DLocalFiniteElement<double,double,4> pk42d;
  testPk(pk42d);

  Pk3DLocalFiniteElement<double,double,1> pk13d;
  testPk(pk13d);
  Pk3DLocalFiniteElement<double,double,2> pk23d;
  testPk(pk23d);
  Pk3DLocalFiniteElement<double,double,3> pk33d;
  testPk(pk33d);
  Pk3DLocalFiniteElement<double,double,4> pk43d;
  testPk(pk43d);

  return success ? 0 : 1;
}
catch (Exception e) {

  std::cerr << e << std::endl;
  return 1;
}
