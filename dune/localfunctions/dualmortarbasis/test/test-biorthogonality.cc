// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <algorithm>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/lagrange/p1.hh>
#include <dune/localfunctions/lagrange/q1.hh>

#include <dune/localfunctions/dualmortarbasis/dualp1.hh>
#include <dune/localfunctions/dualmortarbasis/dualq1.hh>

double TOL = 1e-9;

// \brief Test if dual and Lagrange functions are bi-orthogonal, i.e. \integral dual_i * lagrange_j = \delta_ij \integral lagrange_j
template<class DualLfe,class LagrangeLfe>
bool testBiorthogonality(const DualLfe& dualLfe, const LagrangeLfe& lagrangeLfe)
{
  const int dim = DualLfe::Traits::LocalBasisType::Traits::dimDomain;
  const Dune::QuadratureRule<double,dim>& quad = Dune::QuadratureRules<double,dim>::rule(dualLfe.type(),dualLfe.localBasis().order() + lagrangeLfe.localBasis().order());

  const unsigned int numLagBasFct = lagrangeLfe.localBasis().size();
  const unsigned int numDualBasFct = dualLfe.localBasis().size();

  // save the integrals of all mixed combinations in mixedMassMat[numLagBasFct][numDualBasFct]
  std::vector<std::vector<double> > mixedMassMat(numLagBasFct);
  for (unsigned int i = 0; i < numLagBasFct; i++)
    mixedMassMat[i].resize(numDualBasFct);

  // integrate all lagrange basis functions
  std::vector<double> integralLagrange(numLagBasFct);

  for (unsigned int k=0; k<numLagBasFct; k++)
    std::fill(mixedMassMat[k].begin(), mixedMassMat[k].end(), 0.0);

  std::fill(integralLagrange.begin(), integralLagrange.end(), 0.0);

  for(size_t i=0; i<quad.size(); i++) {
    const Dune::FieldVector<double,dim>& pos = quad[i].position();
    std::vector<typename DualLfe::Traits::LocalBasisType::Traits::RangeType> dualValues(dualLfe.localBasis().size());
    std::vector<typename LagrangeLfe::Traits::LocalBasisType::Traits::RangeType> lagrangeValues(lagrangeLfe.localBasis().size());

    // evaluate basis functions
    dualLfe.localBasis().evaluateFunction(pos,dualValues);
    lagrangeLfe.localBasis().evaluateFunction(pos,lagrangeValues);

    // check if the values sum up to one
    double unity = 0;
    for (unsigned int l=0; l<numDualBasFct; l++)
      unity += dualValues[l];
    if (fabs(unity -1) > TOL)
      std::cout<<"Dual basis functions don't sum up to 1\n";

    double weight = quad[i].weight();

    for (unsigned int k=0; k<numLagBasFct; k++) {

      integralLagrange[k] += lagrangeValues[k]*weight;

      for (unsigned int l=0; l<numDualBasFct; l++)
        mixedMassMat[k][l] += weight*(lagrangeValues[k]*dualValues[l]);

    }
  }
  bool biorthog = true;
  // check if diagonal of mixedMassMat equals the lagrange integral and the off-diagonals are zero
  for (unsigned int k=0; k<numLagBasFct; k++)
    for (unsigned int l=0; l<numLagBasFct; l++)
      if(std::fabs(mixedMassMat[k][l] - (k==l)*(integralLagrange[k])) > TOL) {
        std::cout<<"Dual basis function "<<l<<" and Lagrange basis function "<<k
                 <<" are note bi-orthogonal. The values differs by "
                 <<std::fabs(mixedMassMat[k][l]-(k==l)*integralLagrange[k])<<std::endl;
        biorthog = false;
      }

  return biorthog;
}


int main(int argc, char** argv) try
{
  bool success = true;
  Dune::DualP1LocalFiniteElement<double,double,1> dualP1lfem1D;
  Dune::P1LocalFiniteElement<double,double,1> p1lfem1D;
  success = testBiorthogonality(dualP1lfem1D,p1lfem1D) and success;

  Dune::DualP1LocalFiniteElement<double,double,2> dualP1lfem2D;
  Dune::P1LocalFiniteElement<double,double,2> p1lfem2D;
  success = testBiorthogonality(dualP1lfem2D,p1lfem2D) and success;

  Dune::DualP1LocalFiniteElement<double,double,3> dualP1lfem3D;
  Dune::P1LocalFiniteElement<double,double,3> p1lfem3D;
  success = testBiorthogonality(dualP1lfem3D,p1lfem3D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,1> dualQ1lfem1D;
  Dune::Q1LocalFiniteElement<double,double,1> q1lfem1D;
  success = testBiorthogonality(dualQ1lfem1D,q1lfem1D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,2> dualQ1lfem2D;
  Dune::Q1LocalFiniteElement<double,double,2> q1lfem2D;
  success = testBiorthogonality(dualQ1lfem2D,q1lfem2D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,3> dualQ1lfem3D;
  Dune::Q1LocalFiniteElement<double,double,3> q1lfem3D;
  success = testBiorthogonality(dualQ1lfem3D,q1lfem3D) and success;

  return success ? 0 : 1;
}
catch(Dune::Exception e)
{
  std::cout<<e<<std::endl;
  return 1;
}
