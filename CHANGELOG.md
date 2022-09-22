<!--
SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
-->

# Master (will become release 2.9)

## Deprecations and removals

- Deprecated many of the Lagrange headers, use
  `lagrange(cube|prism|pyramid|simplex).hh` instead.

# Release 2.8

* Passing functions that support `f.evaluate(x,y)` to `interpolate()`
  is deprecated. Instead the functions should now provide `operator()`.
  Passing functions providing the old interface is still supported in 2.8.
  * `LocalFiniteElementFunctionBase` is deprecated. You can rely
    on duck-typing when passing functions with the new interface.
  * The virtual interface for interpolating functions in `LocalFiniteElementVirtualInterface`
    now uses `std::function` instead of the deprecated `VirtualFunction`
    for the passed function.
  * The virtual interface wrapper `LocalFiniteElementVirtualImp` now
    requires that the wrapped `LocalFiniteElement` implementation
    supports the new `operator()` based interpolation interface.

* Add an implementation of the Nédélec element of the first kind,
  as introduced in "Nédélec, Mixed finite elements in R^3, 1980,
  DOI: http://dx.doi.org/10.1007/BF01396415".
  Only the first-order case for triangles, tetrahedra, squares and cubes is implemented.

* Fix a bug in a shape function of the second-order Lagrange element
  on the three-dimensional pyramid.

* Add an implementation of the Raviart-Thomas element for tetrehedra with order 0.

* Remove deprecated `GenericLocalFiniteElement::topologyId()`, use
  `type().id()` instead.

* Imported the Python bindings from the 2.7 branch of dune-python.

* Replaced the combination of function arguments `topologyId` and `dim` with a single `GeometryType` argument.
  Tagged the old versions of: `numLagrangePoints`, `equidistantLagrangePoints`, `RTL2InterpolationBuilder::topologyId()`,
  `VirtualMonomialBasis(topologyId)`, `VirtualMonomialBasis::topologyId()` as deprecated.

* Add a construction algorithm for high order Nédélec elements on triangles and tetrahedra.

# Release 2.7

* The header `lagrange.hh` now includes all headers of all Lagrange implementations,
  not just the ones with run-time order.

* Introduce a run-time polymorphic container `LocalFiniteElementVariant`.
  Much like `std::variant`, it implements a type-safe
  union of different `LocalFiniteElement` implementations.  Elements of type
  `LocalFiniteElementVariant` can hold one object from a list of types
  given as template parameters.  These types must be implementations of
  the `LocalFiniteElement` interface, and the container will in turn
  implement this interface.

  Such a `variant`-based polymorphism is not as flexible as full type erasure,
  but it is much easier to implement.  What is more, it is believed that in
  many situations the restriction to a fixed set of implementation types
  is not a problem.

* Add support for `operator()` syntax to `interpolate()`. All `interpolate()`
  implementations now support functions `f` that either support `f.evaluate(x,y)`
  or `y = f(x)`.

* Add an implementation of the Crouzeix-Raviart element.

* Add an implementation of the Brezzi-Douglas-Fortin-Marini element.
  The coefficients and interpolation are implemented for arbitrary
  dimension (>1) and order (>0). The actual basis is only implemented
  for dim=2 and order=1,2,3.

  See core/dune-localfunctions!105 and core/dune-localfunctions!145

* Introduce a convenience header `hierarchical.hh` that includes
  all hierarchical FE implementations.

* Introduce a new class `LagrangeSimplexLocalFiniteElement`, which implements
  Lagrange finite elements on simplices with compile-time dimension and order.
  It currently does not cover more general dimension/order situations than
  what is already available in dune-localfunctions, but it gathers the
  plethora of different Pk3DNodal, PkNodal, P1Nodal, etc implementations
  under one single name.

* Introduce new class `BrezziDouglasMariniSimplexLocalFiniteElement`
  (and the same for cubes) that subsumes all currently existing simplex
  BDM element implementations under a single name.  Domain dimension and
  polynomial order are template parameters now.

* Introduce a convenience header `dune/localfunctions/brezzidouglasmarini.hh`
  that includes all BDM implementations.

# Release 2.6

*  The `diffOrder` value has disappeared from the `LocalBasisTraits` class.
   This value encoded the highest partial derivative order implemented by
   a local basis. Encoding this value as a compile-time parameter led to
   various problems related to the dynamic interface, mainly because it
   became part of the type of the local finite element.  At the same time,
   it was suspected that very few people ever actually used the parameter.

    More practically, two things have disappeared: the `diffOrder` member
    of the `LocalBasisTraits` class, and the 8th template parameter `dorder`
    of that class.  There is no replacement, and if you have used `diffOrder`
    then you currently have to find a way to live without it.  As mentioned
    we believe that this concerns only a very small number of people.

    If you do use `diffOrder` and you absolutely need it or something similar,
    then we'd like to hear from you.  One of the reasons why there is no
    replacement is that we couldn't really think of a good use case to begin with.

*  The `QkLocalFiniteElement` class implements second partial derivatives
   of shape functions now.

* The `clone()` method was removed from the raw (non-virtual) `LocalFiniteElement`
  implementations. If you want to copy a `LocalFiniteElement` in a portable
  way which works for raw implementations as well as for the virtual interface
  class, you have to replace `lfe.clone()` by
  `Dune::LocalFiniteElementCloneFactory<LFEType>::clone(lfe)`.
