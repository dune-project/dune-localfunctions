// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>

#include <dune/common/gmpfield.hh>

// Lagrange type elements
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

// DG type elements
#include <dune/localfunctions/orthonormal.hh>

// Raviart Thomas type elements
#include <dune/localfunctions/raviartthomas.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/grid/common/mcmgmapper.hh>

const bool verbose = false;
const double TOL = 1e-14;

template<class FE,class GV>
class Func :
  public Dune::LocalFiniteElementFunctionBase<FE>::FunctionBase
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;
  typedef typename GV::template Codim<0>::Entity Element;
  typedef typename Element::Geometry Geometry;

  void evaluate (const DomainType& x, RangeType& y) const
  {
    auto gl = element_.geometry().global(x);
    y = 0;
    y[0] = std::abs(std::sin(17.*gl.two_norm2()))+1.;
  }
  void bind(const Element &element)
  {
    element_ = element;
  }
  void unbind()
  {
  }
private:
  Element element_;
};
// store data on all entities excepts elements and vertices since these
// induce no twist issue
template<int dimgrid> struct AllLayout {
  bool contains (Dune::GeometryType gt)
  { return gt.dim()>0 && gt.dim()<dimgrid; }
};
template <class LFE, class GV>
std::size_t testTwist(GV &gridView, int order)
{
  Func<LFE,GV> f;
  std::size_t failures = 0;
  typedef typename GV::Grid GridType;

  Dune::FieldVector<unsigned int,GridType::dimension> maxDofsPerCodim(0);
  unsigned int maxDofs = 0;
  for (auto element : Dune::elements(gridView))
  {
    LFE lfe(element.type(), order);
    auto coefficients = lfe.localCoefficients();
    for (unsigned int c=1; c<element.dimension;++c)
    {
      for (unsigned int s=0; s<element.subEntities(c); ++s)
      {
        unsigned int numDofs = 0;
        for (unsigned int k=0; k<coefficients.size();++k)
        {
          auto key = coefficients.localKey(k);
          if (key.codim() == c && key.subEntity() == s)
            ++numDofs;
        }
        maxDofsPerCodim[c] = std::max(maxDofsPerCodim[c],numDofs);
      }
      maxDofs = std::max(maxDofs,maxDofsPerCodim[c]);
    }
  }

  Dune::MultipleCodimMultipleGeomTypeMapper<GV,AllLayout> mapper(gridView);
  // value=-1: unused (note that global function is >=1 everywhere)
  std::vector<double> globalValues(mapper.size()*maxDofs, -1);

  for (auto element : Dune::elements(gridView))
  {
    if (verbose)
      std::cout << "new element" << std::endl;
    // interpolate the grid function
    f.bind(element);
    LFE lfe(element.type(), order);

    std::vector<double> localValues;
    lfe.localInterpolation().interpolate(f, localValues);

    auto coefficients = lfe.localCoefficients();
    for (unsigned int c=1; c<element.dimension;++c)
    {
      for (unsigned int s=0; s<element.subEntities(c); ++s)
      {
        // first extract the dofs on the subentity (note: assume always the same)
        std::vector<double> subEntityValues(maxDofsPerCodim[c],-1);
        for (unsigned int k=0; k<coefficients.size();++k)
        {
          auto key = coefficients.localKey(k);
          if (key.codim() == c && key.subEntity() == s)
            subEntityValues[key.index()] = localValues[k];
        }
        assert( subEntityValues.size() <= maxDofs );
        // apply twist matrix here (or possibly inverse)
        // matrix.mv(subEntityValue)
        auto index = maxDofs*mapper.subIndex(element, s, c);
        if (globalValues[index] == -1)
        {
          // first time at this subentity
          for (unsigned int i=0;i<subEntityValues.size();++i)
            globalValues[index+i] = subEntityValues[i];
        }
        else
        {
          // been there so compare values
          std::vector<double> referenceValues(subEntityValues.size());
          for (unsigned int i=0;i<subEntityValues.size();++i)
            referenceValues[i] = globalValues[index+i];
          bool failed = false;
          for (unsigned int i=0;i<subEntityValues.size();++i)
          {
            if (std::abs(referenceValues[i] - subEntityValues[i]) > TOL)
            {
              failed = true;
              ++failures;
              if (verbose)
                std::cout << "    (" << s << "," << c << "): "
                  << "refval[" << i << "]="
                  << referenceValues[i] << " "
                  << "new value[" << i << "]="
                  << subEntityValues[i] << std::endl;
            }
          }
        }
      }
    }
  }
  return failures;
}
template <class GridType>
void testTwist(GridType &grid)
{
  static const int dimension = GridType::dimension;
  for (int level=0; level<2; ++level)
  {
    typedef typename GridType::LeafGridView GridView;
    GridView gridView = grid.leafGridView();
    std::cout << "test on level: " << level << std::endl;
    for (unsigned int order=2; order<=4; ++order)
    {
      std::cout << "order : " << order << std::endl;
      typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,dimension,double,double> LFE;
      std::size_t failures = testTwist<LFE>(gridView,order);
      if (failures > 0)
        std::cout << "detwisting failed" << std::endl;
    }
    grid.globalRefine(1);
  }
}
int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );
  static const int dimension = 3;
  {
    std::cout << "testing YaspGrid - should be twist free" << std::endl;
    typedef Dune::YaspGrid<dimension> GridType;
    Dune::FieldVector<double,dimension> bbox = {1, 1, 1};
    std::array<int,dimension> el = {1, 1, 1};
    GridType grid(bbox,el);
    testTwist(grid);
  }
  {
    std::cout << "testing ALUGrid<cube> using a twist free grid" << std::endl;
    typedef Dune::ALUGrid<dimension,dimension,Dune::cube,Dune::nonconforming> GridType;
    Dune::GridPtr< GridType > gridPtr( "notwist.dgf" );
    GridType &grid = *gridPtr;
    testTwist(grid);
  }
  {
    std::cout << "testing ALUGrid<cube> using a grid with twist" << std::endl;
    typedef Dune::ALUGrid<dimension,dimension,Dune::cube,Dune::nonconforming> GridType;
    Dune::GridPtr< GridType > gridPtr( "twist.dgf" );
    GridType &grid = *gridPtr;
    testTwist(grid);
  }
}
catch (const Dune::Exception &e)
{
  std::cout << e << std::endl;
  return 1;
}
