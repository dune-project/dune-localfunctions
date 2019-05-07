# Master (will become release 2.7)

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

* Introduce a convenience header `hierarchical.hh` that includes
  all hierarchical FE implementations.

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
