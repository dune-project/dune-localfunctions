// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

// This header is not part of the official Dune API and might be subject
// to change.  You can use this header to test external finite element
// implementations, but be warned that your tests might break with future
// Dune versions.

#ifndef DUNE_LOCALFUNCTIONS_TEST_GEOMETRIES_HH
#define DUNE_LOCALFUNCTIONS_TEST_GEOMETRIES_HH

#include <cstddef>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

template<class ctype, std::size_t dim>
class TestGeometries;

template<class ctype>
class TestGeometries<ctype, 0> :
  public std::vector<Dune::MultiLinearGeometry<ctype, 0, 0> >
{
  static const std::size_t dim = 0;

public:
  typedef Dune::MultiLinearGeometry<ctype, dim, dim> Geometry;

  TestGeometries() {
    Dune::GeometryType gt;
    std::vector<Dune::FieldVector<ctype, dim> > coords;

    gt = Dune::GeometryTypes::vertex;
    coords.resize(1);
    this->push_back(Geometry(gt, coords));
  }

  const Geometry &get(const Dune::GeometryType &gt) const {
    for(std::size_t i = 0; i < this->size(); ++i)
      if((*this)[i].type() == gt) return (*this)[i];
    DUNE_THROW(Dune::NotImplemented, "No predefined test-geometry in "
               "dimension " << dim << " for GeometryType " << gt);
  }
};

template<class ctype>
class TestGeometries<ctype, 1> :
  public std::vector<Dune::MultiLinearGeometry<ctype, 1, 1> >
{
  static const std::size_t dim = 1;

public:
  typedef Dune::MultiLinearGeometry<ctype, dim, dim> Geometry;

  TestGeometries() {
    Dune::GeometryType gt;
    std::vector<Dune::FieldVector<ctype, dim> > coords;

    gt = Dune::GeometryTypes::line;
    coords.resize(2);
    coords[0][0] = -.3;
    coords[1][0] =  .7;
    this->push_back(Geometry(gt, coords));
  }

  const Geometry &get(const Dune::GeometryType &gt) const {
    for(std::size_t i = 0; i < this->size(); ++i)
      if((*this)[i].type() == gt) return (*this)[i];
    DUNE_THROW(Dune::NotImplemented, "No predefined test-geometry in "
               "dimension " << dim << " for GeometryType " << gt);
  }
};

template<class ctype>
class TestGeometries<ctype, 2> :
  public std::vector<Dune::MultiLinearGeometry<ctype, 2, 2> >
{
  static const std::size_t dim = 2;

public:
  typedef Dune::MultiLinearGeometry<ctype, dim, dim> Geometry;

  TestGeometries() {
    Dune::GeometryType gt;
    std::vector<Dune::FieldVector<ctype, dim> > coords;

    gt = Dune::GeometryTypes::triangle;
    coords.resize(3);
    coords[0][0] = -.5; coords[0][1] = -.5;
    coords[1][0] =  .5; coords[1][1] = -.5;
    coords[2][0] = 0  ; coords[2][1] =  .5;
    this->push_back(Geometry(gt, coords));

    gt = Dune::GeometryTypes::quadrilateral;
    coords.resize(4);
    coords[0][0] = -.5; coords[0][1] = 0;
    coords[1][0] = 0  ; coords[1][1] = -.5;
    coords[2][0] =  .5; coords[2][1] = 0;
    coords[3][0] = 0  ; coords[3][1] =  .5;
    this->push_back(Geometry(gt, coords));
  }

  const Geometry &get(const Dune::GeometryType &gt) const {
    for(std::size_t i = 0; i < this->size(); ++i)
      if((*this)[i].type() == gt) return (*this)[i];
    DUNE_THROW(Dune::NotImplemented, "No predefined test-geometry in "
               "dimension " << dim << " for GeometryType " << gt);
  }
};

template<class ctype>
class TestGeometries<ctype, 3> :
  public std::vector<Dune::MultiLinearGeometry<ctype, 3, 3> >
{
  static const std::size_t dim = 3;

public:
  typedef Dune::MultiLinearGeometry<ctype, dim, dim> Geometry;

  TestGeometries() {
    Dune::GeometryType gt;
    std::vector<Dune::FieldVector<ctype, dim> > coords;

    gt = Dune::GeometryTypes::tetrahedron;
    coords.resize(4);
    coords[0][0] = -.5; coords[0][1] = -.5; coords[0][2] = -.5;
    coords[1][0] =  .5; coords[1][1] = -.5; coords[1][2] = -.5;
    coords[2][0] = 0  ; coords[2][1] =  .5; coords[2][2] = -.5;
    coords[3][0] = 0  ; coords[3][1] =  0 ; coords[3][2] =  .5;
    this->push_back(Geometry(gt, coords));

    gt = Dune::GeometryTypes::pyramid;
    coords.resize(5);
    coords[0][0] = -.5; coords[0][1] = 0;   coords[0][2] = -.5;
    coords[1][0] = 0  ; coords[1][1] = -.5; coords[1][2] = -.5;
    coords[2][0] =  .5; coords[2][1] = 0;   coords[2][2] = -.5;
    coords[3][0] = 0  ; coords[3][1] =  .5; coords[3][2] = -.5;
    coords[4][0] =  .1; coords[4][1] =  .1; coords[4][2] =  .1;
    this->push_back(Geometry(gt, coords));

    gt = Dune::GeometryTypes::prism;
    coords.resize(6);
    coords[0][0] = -.6; coords[0][1] = -.5; coords[0][2] = -.4;
    coords[1][0] =  .5; coords[1][1] = -.6; coords[1][2] = -.5;
    coords[2][0] =  .1; coords[2][1] =  .5; coords[2][2] = -.6;
    coords[3][0] = -.5; coords[3][1] = -.4; coords[3][2] =  .5;
    coords[4][0] =  .4; coords[4][1] = -.5; coords[4][2] =  .6;
    coords[5][0] = 0  ; coords[5][1] =  .4; coords[5][2] =  .5;
    this->push_back(Geometry(gt, coords));

    gt = Dune::GeometryTypes::hexahedron;
    coords.resize(8);
    coords[0][0] = -.7; coords[0][1] = -.6; coords[0][2] = -.5;
    coords[1][0] =  .4; coords[1][1] = -.3; coords[1][2] = -.7;
    coords[2][0] = -.6; coords[2][1] =  .5; coords[2][2] = -.4;
    coords[3][0] =  .3; coords[3][1] =  .7; coords[3][2] = -.6;
    coords[4][0] = -.5; coords[4][1] = -.4; coords[4][2] =  .3;
    coords[5][0] =  .7; coords[5][1] = -.6; coords[5][2] =  .5;
    coords[6][0] = -.4; coords[6][1] =  .3; coords[6][2] =  .7;
    coords[7][0] =  .6; coords[7][1] =  .5; coords[7][2] =  .4;
    this->push_back(Geometry(gt, coords));
  }

  const Geometry &get(const Dune::GeometryType &gt) const {
    for(std::size_t i = 0; i < this->size(); ++i)
      if((*this)[i].type() == gt) return (*this)[i];
    DUNE_THROW(Dune::NotImplemented, "No predefined test-geometry in "
               "dimension " << dim << " for GeometryType " << gt);
  }
};

#endif // DUNE_LOCALFUNCTIONS_TEST_GEOMETRIES_HH
