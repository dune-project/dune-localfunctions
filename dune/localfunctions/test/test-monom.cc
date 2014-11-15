// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/monomial.hh>

#include "geometries.hh"
#include "test-fe.hh"

// tolerance for floating-point comparisons
const double eps = 1e-9;
// stepsize for numerical differentiation
const double delta = 1e-5;

template<int dim>
struct Dim {
  template<int p>
  struct Order {

    static void apply(int &result) {
      std::cout << "== Checking global-valued Monom elements (with "
                << "dim=" << dim << ", p=" << p << ")" << std::endl;

      typedef TestGeometries<double, dim> TestGeos;
      static const TestGeos testGeos;

      typedef typename TestGeos::Geometry Geometry;
      Dune::MonomFiniteElementFactory<Geometry, double, p> feFactory;

      for(std::size_t i = 0; i < testGeos.size(); ++i) {
        const Geometry &geo = testGeos[i];

        std::cout << "=== GeometryType " << geo.type() << std::endl;

        bool success = testFE(geo, feFactory.make(geo), eps, delta);

        if(success && result != 1)
          result = 0;
        else
          result = 1;
      }
    }
  };

  static void apply(int &result)
  { Dune::ForLoop<Order, 0, 3>::apply(result); }
};

int main(int argc, char** argv) {
  try {
    int result = 77;

    Dune::ForLoop<Dim, 1, 3>::apply(result);

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
