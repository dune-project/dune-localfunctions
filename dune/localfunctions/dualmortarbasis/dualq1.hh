// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DUAL_Q1_LOCALFINITEELEMENT_HH
#define DUNE_DUAL_Q1_LOCALFINITEELEMENT_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "dualq1/dualq1localbasis.hh"
#include "dualq1/dualq1localcoefficients.hh"
#include "dualq1/dualq1localinterpolation.hh"

namespace Dune
{

  /** \brief The local dual Q1 finite element on cubes
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the hypercube
   */
  template<class D, class R, int dim>
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
      gt.makeCube(dim);

      // dual basis functions are linear combinations of Lagrange elements
      // compute these coefficients here because the basis and the local interpolation needs them

      const Dune::QuadratureRule<D,dim>& quad = Dune::QuadratureRules<D,dim>::rule(gt, 2*dim);

      const int size = 1 <<dim;

      // assemble mass matrix on the reference element
      Dune::FieldMatrix<R, size, size> massMat;
      massMat = 0;

      // and the integrals of the lagrange shape functions
      std::vector<Dune::FieldVector<R,1> > integral(size);
      for (int i=0; i<size; i++)
        integral[i] = 0;

      for(size_t pt=0; pt<quad.size(); pt++) {

        const Dune::FieldVector<D ,dim>& pos = quad[pt].position();
        std::vector<typename Traits::LocalBasisType::Traits::RangeType> q1Values(size);

        // evaluate q1 basis functions
        for (int i=0; i<size; i++) {

          q1Values[i] = 1;

          for (int j=0; j<dim; j++)
            // if j-th bit of i is set multiply with in[j], else with 1-in[j]
            q1Values[i] *= (i & (1<<j)) ? pos[j] :  1-pos[j];
        }

        double weight = quad[pt].weight();

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
      Dune::array<Dune::FieldVector<R, size>, size> coefficients;

      for (int i=0; i<size; i++) {

        Dune::FieldVector<R, size> rhs(0);
        rhs[i] = integral[i];

        coefficients[i] = 0;
        massMat.solve(coefficients[i] ,rhs);

      }

      basis.setCoefficients(coefficients);
      interpolation.setCoefficients(coefficients);
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

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    DualQ1LocalFiniteElement* clone () const
    {
      return new DualQ1LocalFiniteElement(*this);
    }

  private:
    DualQ1LocalBasis<D,R,dim> basis;
    DualQ1LocalCoefficients<dim> coefficients;
    DualQ1LocalInterpolation<dim,DualQ1LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };
}

#endif
