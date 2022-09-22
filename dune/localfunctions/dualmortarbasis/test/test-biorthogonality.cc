// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <algorithm>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include <dune/localfunctions/dualmortarbasis/dualp1.hh>
#include <dune/localfunctions/dualmortarbasis/dualq1.hh>

double TOL = 1e-9;

//! \brief Test if dual and Lagrange functions are bi-orthogonal, i.e. \integral dual_i * lagrange_j = \delta_ij \integral lagrange_j
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


//! \brief Test if dual and Lagrange functions are bi-orthogonal on adjacent faces i.e. \integral dual_i * lagrange_j = \delta_ij \integral lagrange_j for each face that i is part of.
template<class DualLfe,class LagrangeLfe>
bool testFaceBiorthogonality(const DualLfe& dualLfe, const LagrangeLfe& lagrangeLfe)
{
  using DualTraits = typename DualLfe::Traits::LocalBasisType::Traits;
  using LagrangeTraits = typename LagrangeLfe::Traits::LocalBasisType::Traits;

  const int dim = DualTraits::dimDomain;
  assert (dim == LagrangeTraits::dimDomain);

  using ctype = typename DualTraits::DomainFieldType;
  using field_type = typename DualTraits::RangeFieldType;

  const auto& refElement = Dune::ReferenceElements<ctype,dim>::general(dualLfe.type());

  const unsigned int numLagBasFct = lagrangeLfe.localBasis().size();
  const unsigned int numDualBasFct = dualLfe.localBasis().size();

  // save the integrals of all mixed combinations
  std::vector<std::vector<field_type> > mixedMassMat(numLagBasFct);
  for (unsigned int i = 0; i < numLagBasFct; i++)
      mixedMassMat[i].resize(numDualBasFct);

  // integrate all lagrange basis functions
  std::vector<field_type> integralLagrange(numLagBasFct);

  std::vector<typename DualTraits::RangeType> dualValues(dualLfe.localBasis().size());
  std::vector<typename LagrangeTraits::RangeType> lagrangeValues(lagrangeLfe.localBasis().size());

  bool biorthog = true;
  // loop over the number of faces
  for (int i=0; i<refElement.size(1); i++) {

    // get geometry
    const auto& geometry = refElement.template geometry<1>(i);
    // get a quadrature rule
    const auto& quad = Dune::QuadratureRules<ctype,dim-1>::rule(refElement.type(i,1),  dualLfe.localBasis().order() + lagrangeLfe.localBasis().order());

    // initialise with zeroes
    for (unsigned int k=0; k<numLagBasFct; k++)
      std::fill(mixedMassMat[k].begin(), mixedMassMat[k].end(), 0.0);

    std::fill(integralLagrange.begin(), integralLagrange.end(), 0.0);

    for(size_t j=0; j < quad.size(); j++) {

      const auto& pos = quad[j].position();
      const auto& elementPos = geometry.global(pos);

      // evaluate basis functions
      dualLfe.localBasis().evaluateFunction(elementPos,dualValues);
      lagrangeLfe.localBasis().evaluateFunction(elementPos,lagrangeValues);

      auto weight = quad[j].weight();

      for (unsigned int k=0; k<numLagBasFct; k++) {

        integralLagrange[k] += lagrangeValues[k]*weight;

        for (unsigned int l=0; l<numDualBasFct; l++)
          mixedMassMat[k][l] += weight*(lagrangeValues[k]*dualValues[l]);
      }
    }

    // check if diagonal of mixedMassMat equals the lagrange integral and the off-diagonals are zero
    for (unsigned int p=0; p<(unsigned int)refElement.size(i,1,dim); p++) {
      int k= refElement.subEntity(i,1,p,dim);
      for (unsigned int q=0; q<(unsigned int)refElement.size(i,1,dim); q++) {
        int l= refElement.subEntity(i,1,q,dim);
        if(std::fabs(mixedMassMat[k][l] - (k==l)*(integralLagrange[k])) > TOL) {
          std::cout<<"for face "<<i<<" : \n";
          std::cout<<"Dual basis function "<<l<<" and Lagrange basis function "<<k
            <<" are not bi-orthogonal. The values differs by "
            <<std::fabs(mixedMassMat[k][l]-(k==l)*integralLagrange[k])<<std::endl;
          biorthog = false;
        }
      }
    }
  }
  return biorthog;
}

int main(int argc, char** argv)
{
  bool success = true;
  Dune::DualP1LocalFiniteElement<double,double,1> dualP1lfem1D;
  Dune::LagrangeSimplexLocalFiniteElement<double,double,1,1> p1lfem1D;
  success = testBiorthogonality(dualP1lfem1D,p1lfem1D) and success;

  Dune::DualP1LocalFiniteElement<double,double,2> dualP1lfem2D;
  Dune::LagrangeSimplexLocalFiniteElement<double,double,2,1> p1lfem2D;
  success = testBiorthogonality(dualP1lfem2D,p1lfem2D) and success;

  Dune::DualP1LocalFiniteElement<double,double,3> dualP1lfem3D;
  Dune::LagrangeSimplexLocalFiniteElement<double,double,3,1> p1lfem3D;
  success = testBiorthogonality(dualP1lfem3D,p1lfem3D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,1> dualQ1lfem1D;
  Dune::LagrangeCubeLocalFiniteElement<double,double,1,1> q1lfem1D;
  success = testBiorthogonality(dualQ1lfem1D,q1lfem1D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,2> dualQ1lfem2D;
  Dune::LagrangeCubeLocalFiniteElement<double,double,2,1> q1lfem2D;
  success = testBiorthogonality(dualQ1lfem2D,q1lfem2D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,3> dualQ1lfem3D;
  Dune::LagrangeCubeLocalFiniteElement<double,double,3,1> q1lfem3D;
  success = testBiorthogonality(dualQ1lfem3D,q1lfem3D) and success;

  Dune::DualP1LocalFiniteElement<double,double,1,true> dualFaceP1lfem1D;
  success = testFaceBiorthogonality(dualFaceP1lfem1D,p1lfem1D) and success;

  Dune::DualP1LocalFiniteElement<double,double,2,true> dualFaceP1lfem2D;
  success = testFaceBiorthogonality(dualFaceP1lfem2D,p1lfem2D) and success;

  Dune::DualP1LocalFiniteElement<double,double,3,true> dualFaceP1lfem3D;
  success = testFaceBiorthogonality(dualFaceP1lfem3D,p1lfem3D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,1,true> dualFaceQ1lfem1D;
  success = testFaceBiorthogonality(dualFaceQ1lfem1D,q1lfem1D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,2,true> dualFaceQ1lfem2D;
  success = testFaceBiorthogonality(dualFaceQ1lfem2D,q1lfem2D) and success;

  Dune::DualQ1LocalFiniteElement<double,double,3,true> dualFaceQ1lfem3D;
  success = testFaceBiorthogonality(dualFaceQ1lfem3D,q1lfem3D) and success;

  return success ? 0 : 1;
}
