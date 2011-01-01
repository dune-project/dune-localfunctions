// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/utility/mockgeometry.hh>
#include <dune/grid/utility/vertexorder.hh>

#include <dune/localfunctions/whitney/edges0.5.hh>

#include "test-fe.hh"

int main(int argc, char** argv) {
  try {
    int result = 77;

    // tolerance for floating-point comparisons
    static const double eps = 1e-9;
    // stepsize for numerical differentiation
    static const double delta = 1e-5;

    { // 2D
      std::cout << "== Checking global-valued EdgeS0_5 elements (with dim=2)"
                << std::endl;

      Dune::GeometryType gt;
      gt.makeTriangle();

      Dune::FieldVector<double, 2> corners[3];
      corners[0][0] = -.5; corners[0][1] = -.5;
      corners[1][0] =  .5; corners[1][1] = -.5;
      corners[2][0] = 0  ; corners[2][1] =  .5;
      typedef Dune::MockGeometry<double, 2, 2> Geometry;
      Geometry geo(gt, corners);

      std::size_t vertexIds[] = {0, 1, 2};
      Dune::GeneralVertexOrder<2, std::size_t>
      vo(gt, vertexIds+0, vertexIds+3);

      Dune::EdgeS0_5FiniteElementFactory<Geometry, double> feFactory;
      bool success = testFE(geo, feFactory.make(geo, vo), eps, delta);

      if(success && result != 1)
        result = 0;
      else
        result = 1;
    }

    { // 3D
      std::cout << "== Checking global-valued EdgeS0_5 elements (with dim=3)"
                << std::endl;

      Dune::GeometryType gt;
      gt.makeTetrahedron();

      Dune::FieldVector<double, 3> corners[4];
      corners[0][0] = -.5; corners[0][1] = -.5; corners[0][2] = -.5;
      corners[1][0] =  .5; corners[1][1] = -.5; corners[1][2] = -.5;
      corners[2][0] = 0  ; corners[2][1] =  .5; corners[2][2] = -.5;
      corners[3][0] = 0  ; corners[3][1] =  0 ; corners[3][2] =  .5;
      typedef Dune::MockGeometry<double, 3, 3> Geometry;
      Geometry geo(gt, corners);

      std::size_t vertexIds[] = {0, 1, 2, 3};
      Dune::GeneralVertexOrder<3, std::size_t>
      vo(gt, vertexIds+0, vertexIds+4);

      Dune::EdgeS0_5FiniteElementFactory<Geometry, double> feFactory;
      bool success = testFE(geo, feFactory.make(geo, vo), eps, delta);

      if(success && result != 1)
        result = 0;
      else
        result = 1;
    }

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    throw;
  }
}
