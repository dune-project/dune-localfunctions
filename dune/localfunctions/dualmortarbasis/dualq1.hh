// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_Q1_LOCALFINITEELEMENT_HH
#define DUNE_DUAL_Q1_LOCALFINITEELEMENT_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include "dualq1/dualq1localbasis.hh"
#include "dualq1/dualq1localcoefficients.hh"
#include "dualq1/dualq1localinterpolation.hh"

namespace Dune
{
  /**
   * \brief The local dual Q1 finite element on cubes
   *
   *    Note that if the dual functions are chosen to be dual on the faces,
   *    the integrated product of a Lagrange \f$\lambda_p\f$ and dual
   *    function \f$\theta_q\f$ over faces not containing \f$q\f$ does in
   *    general not vanish.
   *
   * \ingroup DualMortar
   *
   * \tparam D Domain data type
   * \tparam R Range data type
   * \tparam dim Dimension of the hypercube
   * \tparam faceDual If set, the basis functions are bi-orthogonal only on faces containing the corresponding vertex.
   */
  template<class D, class R, int dim, bool faceDual=false>
  class DualQ1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<DualQ1LocalBasis<D,R,dim>,DualQ1LocalCoefficients<dim>,
        DualQ1LocalInterpolation<dim,DualQ1LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    DualQ1LocalFiniteElement ()
    {
      if (faceDual)
          setupFaceDualCoefficients();
      else
          setupDualCoefficients();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    /** \brief Setup the coefficients for the face dual basis functions. */
    void setupFaceDualCoefficients();

    /** \brief Setup the coefficients for the dual basis functions. */
    void setupDualCoefficients();

    DualQ1LocalBasis<D,R,dim> basis;
    DualQ1LocalCoefficients<dim> coefficients;
    DualQ1LocalInterpolation<dim,DualQ1LocalBasis<D,R,dim> > interpolation;
  };

  template<class D, class R, int dim, bool faceDual>
  void DualQ1LocalFiniteElement<D,R,dim,faceDual>::setupDualCoefficients()
  {

    const int size = 1 <<dim;
    std::array<Dune::FieldVector<R, size>, size> coeffs;

    // dual basis functions are linear combinations of Lagrange elements
    // compute these coefficients here because the basis and the local interpolation needs them
    const auto& quad = Dune::QuadratureRules<D,dim>::rule(type(), 2*dim);

    // assemble mass matrix on the reference element
    Dune::FieldMatrix<R, size, size> massMat;
    massMat = 0;

    // and the integrals of the lagrange shape functions
    std::vector<Dune::FieldVector<R,1> > integral(size);
    for (int i=0; i<size; i++)
      integral[i] = 0;

    Dune::Impl::LagrangeCubeLocalBasis<D,R,dim,1> q1Basis;
    for(size_t pt=0; pt<quad.size(); pt++) {

      const Dune::FieldVector<D ,dim>& pos = quad[pt].position();
      std::vector<Dune::FieldVector<R,1> > q1Values(size);
      q1Basis.evaluateFunction(pos,q1Values);

      D weight = quad[pt].weight();

      for (int k=0; k<size; k++) {
        integral[k] += q1Values[k]*weight;

        for (int l=0; l<=k; l++)
          massMat[k][l] += weight*(q1Values[k]*q1Values[l]);
      }
    }

    // make matrix symmetric
    for (int i=0; i<size-1; i++)
      for (int j=i+1; j<size; j++)
        massMat[i][j] = massMat[j][i];

    //solve for the coefficients

    for (int i=0; i<size; i++) {

      Dune::FieldVector<R, size> rhs(0);
      rhs[i] = integral[i];

      coeffs[i] = 0;
      massMat.solve(coeffs[i] ,rhs);

    }

    basis.setCoefficients(coeffs);
    interpolation.setCoefficients(coeffs);
  }

  template<class D, class R, int dim, bool faceDual>
  void DualQ1LocalFiniteElement<D,R,dim,faceDual>::setupFaceDualCoefficients()
  {

    const int size = 1 <<dim;
    std::array<Dune::FieldVector<R, size>, size> coeffs;

    // dual basis functions are linear combinations of Lagrange elements
    Dune::Impl::LagrangeCubeLocalBasis<D,R,dim,1> q1Basis;

    const auto& refElement = Dune::ReferenceElements<D,dim>::general(type());

    // loop over faces
    for (int i=0; i<refElement.size(1);i++) {

      const auto& quad = Dune::QuadratureRules<D,dim-1>::rule(refElement.type(i,1),2*dim);

      // for each face assemble the mass matrix over the face of all
      // non-vanishing basis functions,
      // for cubes refElement.size(i,1,dim)=size/2
      Dune::FieldMatrix<R, size/2, size/2> massMat;
      massMat = 0;

      // get geometry
      const auto& geometry = refElement.template geometry<1>(i);

      // and the integrals of the lagrange shape functions
      std::vector<Dune::FieldVector<R,1> > integral(size/2);
      for (int k=0; k<size/2; k++)
        integral[k] = 0;

      for(size_t pt=0; pt<quad.size(); pt++) {

        const auto& pos = quad[pt].position();
        const auto& elementPos = geometry.global(pos);

        std::vector<Dune::FieldVector<R,1> > q1Values(size);
        q1Basis.evaluateFunction(elementPos,q1Values);

        D weight = quad[pt].weight();

        for (int k=0; k<refElement.size(i,1,dim); k++) {
          int row = refElement.subEntity(i,1,k,dim);
          integral[k] += q1Values[row]*weight;

          for (int l=0; l<refElement.size(i,1,dim); l++) {
            int col = refElement.subEntity(i,1,l,dim);
            massMat[k][l] += weight*(q1Values[row]*q1Values[col]);
          }
        }
      }

      // solve for the coefficients
      // note that we possibly overwrite coefficients for neighbouring faces
      // this is okay since the coefficients are symmetric
      for (int l=0; l<refElement.size(i,1,dim); l++) {

        int row = refElement.subEntity(i,1,l,dim);
        Dune::FieldVector<R, size/2> rhs(0);
        rhs[l] = integral[l];

        Dune::FieldVector<R, size/2> x(0);
        massMat.solve(x ,rhs);

        for (int k=0; k<refElement.size(i,1,dim); k++) {
          int col = refElement.subEntity(i,1,k,dim);
          coeffs[row][col]=x[k];
        }
      }
    }

    basis.setCoefficients(coeffs);
    interpolation.setCoefficients(coeffs);
  }
}
#endif
