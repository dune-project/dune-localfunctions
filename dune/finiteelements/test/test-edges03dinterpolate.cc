// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/geometry.hh>

#include "../edges03d.hh"

class U
{
public:
  template<typename DF, typename RF>
  void evaluate(const Dune::FieldVector<DF, 3> &in,
                Dune::FieldVector<RF, 3> &out) const
  {
    out[0] = 1;
    out[1] = 0;
    out[2] = 0;
  }
};


template<typename ctype>
class UnitTetrahedronGeometry
  : public Dune::GenericGeometry::BasicGeometry<
        3,
        Dune::GenericGeometry::DefaultGeometryTraits<
            ctype,
            3,
            3,
            true
            >
        >
{
  typedef Dune::GenericGeometry::BasicGeometry<
      3,
      Dune::GenericGeometry::DefaultGeometryTraits<
          ctype,
          3,
          3,
          true
          >
      > Base;

  static Dune::FieldVector<Dune::FieldVector<double,3>,4>
  getCoords()
  {
    Dune::FieldVector<Dune::FieldVector<double,3>,4> coords;
    coords[0] = 0;
    coords[1] = 0; coords[1][0] = 1;
    coords[2] = 0; coords[2][1] = 1;
    coords[3] = 0; coords[3][2] = 1;
    return coords;
  }

public:
  UnitTetrahedronGeometry()
    : Base(Dune::GeometryType(Dune::GeometryType::simplex, 3), getCoords())
  {}
};

int main (int argc, char *argv[])
try
{
  typedef Dune::EdgeS03DLocalFiniteElement<double, double> LFE;
  // default orientation
  LFE lfe;

  U u;

  typedef std::vector<double> C;
  C expected(6);
  expected[0] = 1;
  expected[1] = 0;
  expected[2] = -std::sqrt(0.5);
  expected[3] = 0;
  expected[4] = -std::sqrt(0.5);
  expected[5] = 0;

  C interpolated;

  lfe.localInterpolation().interpolateGlobal(u, interpolated,
                                             UnitTetrahedronGeometry<double>());

  for(unsigned i = 0; i < expected.size(); ++i)
    std::cout << "expected[" << i << "] = " << expected[i] << std::endl;

  for(unsigned i = 0; i < interpolated.size(); ++i)
    std::cout << "interpolated[" << i << "] = " << interpolated[i] << std::endl;

  if(expected.size() != interpolated.size()) {
    std::cerr << "Error: Size of expected and interpolated vectors don't match!" << std::endl;
    return 1;
  }

  int result = 0;

  for(unsigned i = 0; i < expected.size(); ++i)
    if(Dune::FloatCmp::ne(expected[i], interpolated[i])) {
      std::cerr << "Error: DoF #" << i << " does not match expected value!" << std::endl;
      result = 1;
    }

  return result;
}
catch (Dune::Exception e) {

  std::cout << e << std::endl;
  return 1;
}
