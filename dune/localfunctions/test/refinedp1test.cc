// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#if HAVE_ALUGRID
#include <cstddef>
#include <iostream>

#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/refined/refinedp1.hh>

template<int dim>
int testForDim()
{
  // //////////////////////////////////////////////////////
  // Generate a grid consisting of a single simplex
  // //////////////////////////////////////////////////////

  typedef typename Dune::ALUGrid<dim,dim,Dune::ALUGridElementType::simplex,
      Dune::ALUGridRefinementType::nonconforming> GridType;

  Dune::GridFactory<GridType> gridFactory;

  Dune::FieldVector<double,dim> pos(0);
  gridFactory.insertVertex(pos);

  for (int i=0; i<dim; i++) {
    pos = 0;   pos[i] = 1;
    gridFactory.insertVertex(pos);
  }

  std::vector<unsigned int> element(dim+1);
  for (int i=0; i<dim+1; i++)
    element[i] = i;

  gridFactory.insertElement( Dune::GeometryType(Dune::GeometryType::simplex,dim), element);

  GridType* gridptr = gridFactory.createGrid();

  GridType& grid = *gridptr;

  // //////////////////////////////////////////////////////
  //   The actual test
  // //////////////////////////////////////////////////////

  typedef typename GridType::template Codim<0>::LeafIterator ElementIterator;
  typedef typename GridType::template Codim<dim>::LeafIterator NodeIterator;

  for (int i=0; i<1; ++i)
    grid.globalRefine(1);

  typename Dune::template PkLocalFiniteElement<typename GridType::ctype,double,dim,1> p1element;
  typename Dune::template RefinedP1LocalFiniteElement<typename GridType::ctype,double,dim> refp1element;
  typename Dune::template PkLocalFiniteElement<typename GridType::ctype,double,dim,2> p2element;

  std::vector<int> indexmap(grid.size(dim),-1);
  std::vector<int> indexmapp2(grid.size(dim),-1);
  std::cout << "test for refinedP1LocalFiniteElement in " << dim << "D ... ";

  for (NodeIterator node=grid.template leafbegin<dim>(); node != grid.template leafend<dim>(); ++node)
  {
    std::vector<typename Dune::FieldVector<double, 1> > refp1values(refp1element.localBasis().size(), 0.0);
    std::vector<typename Dune::FieldVector<double, 1> > p2values(p2element.localBasis().size(), 0.0);
    refp1element.localBasis().evaluateFunction(node->geometry().corner(0),refp1values);
    p2element.localBasis().evaluateFunction(node->geometry().corner(0),p2values);
    int count = 0;
    for (std::size_t i=0; i<refp1element.localBasis().size(); ++i)
    {
      assert(refp1values[i]>=0 and (fabs(1-refp1values[i])<1e-8 or refp1values[i]<1e-8));
      if (fabs(1-refp1values[i])<=1e-8)
      {
        assert(fabs(1-p2values[i])<=1e-8);          // check if ordering of shapefunction is consistent with LocalCoefficients given. As RefinedP1LocalBasis uses the LocalCoefficients of the Pk3DLocalBasis<...2> we compare to the ordering of that basis.
        ++count;
        indexmap[i] = grid.leafIndexSet().template index<dim>(*node);
      }
    }
    assert(count==1);     // check if indexmap is injective meaning one and only one refinedP1LocalBasisFunction evaluates to 1 at current node of refined grid
  }

  for (std::size_t k=0; k<indexmap.size(); ++k)   // check if indexmap is really a map meaning that all refinedP1basisFunctions find a corresponding refined node
  {
    assert(indexmap[k] >= 0);
  }

  for (ElementIterator element = grid.template leafbegin<0>(); element != grid.template leafend<0>(); ++element)
  {
    for (int i=0; i<10; ++i)
    {
      Dune::FieldVector<double,dim> randomPoint;
      do {
        for (int j=0; j<dim; ++j)
          randomPoint[j] = (1.0*rand())/RAND_MAX;
      } while (randomPoint.one_norm() > 1.0);

      std::vector<Dune::FieldVector<double,1> > p1values_,
                                                p1values,
                                                refp1values;

      typedef typename Dune::PkLocalFiniteElement<typename GridType::ctype,double,dim,1>::Traits::LocalBasisType::Traits::JacobianType PkJacobianType;

      std::vector<typename Dune::RefinedP1LocalFiniteElement<typename GridType::ctype,double,dim>::Traits::LocalBasisType::Traits::JacobianType> refp1grads;
      std::vector<PkJacobianType> p1grads_, p1grads;

      p1element.localBasis().evaluateFunction(randomPoint, p1values_);
      refp1element.localBasis().evaluateFunction(element->geometry().global(randomPoint),refp1values);

      p1element.localBasis().evaluateJacobian(randomPoint, p1grads_);
      refp1element.localBasis().evaluateJacobian(element->geometry().global(randomPoint), refp1grads);

      p1values.assign(refp1values.size(),0.0);

      PkJacobianType zeroJacobian;
      zeroJacobian = 0.0;
      p1grads.assign(refp1grads.size(),zeroJacobian);

      const Dune::FieldMatrix<double,dim,dim>& invJacobian = element->geometry().jacobianInverseTransposed(randomPoint);

      for (int v=0; v < element->geometry().corners(); ++v)
      {
        p1values[grid.leafIndexSet().subIndex(*element,v,dim)] = p1values_[v];
        invJacobian.mv(p1grads_[v][0],p1grads[grid.leafIndexSet().subIndex(*element,v,dim)][0]);
      }

      for (std::size_t j=0; j<p1values.size(); ++j)
      {
        assert(fabs(refp1values[j]-p1values[indexmap[j]])<1e-14);         // check if RefinedP1LocalBasis functions on the big tetrahedron evaluate to the same values as the P1LocalBasis functions on refined grid
        for (int k=0; k<dim; ++k)
          assert(fabs(refp1grads[j][0][k]-p1grads[indexmap[j]][0][k])<1e-14);           // do the same for the gradients
      }
    }
  }
  std::cout << "passed." << std::endl;
  return 0;
}


int main (int argc, char** argv) try
{
  //Init MPI
  Dune::MPIHelper::instance(argc, argv);

  testForDim<2>();
  testForDim<3>();
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
}

#endif // HAVE_ALUGRID
