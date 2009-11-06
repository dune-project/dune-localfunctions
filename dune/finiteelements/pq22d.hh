// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQ22DLOCALFINITEELEMENT_HH
#define DUNE_PQ22DLOCALFINITEELEMENT_HH

#include "common/virtualinterface.hh"
#include "q22d.hh"
#include "pk2d.hh"

namespace Dune
{
  template<class D, class R>
  class PQ22DLocalFiniteElement :
    public LocalFiniteElementInterface<
        LocalFiniteElementTraits<
            C1LocalBasisVirtualInterface<
                C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,
                    R,1,Dune::FieldVector<R,1>,
                    Dune::FieldVector<Dune::FieldVector<R,2>,1> > >,
            LocalCoefficientsVirtualInterface,
            LocalInterpolationVirtualInterface< Dune::FieldVector<D,2>, Dune::FieldVector<R,1> >
            >, PQ22DLocalFiniteElement<D,R>
        >
  {
    typedef Dune::FieldVector<D,2> Domain;
    typedef Dune::FieldVector<R,1> Range;
    typedef C1LocalBasisTraits<D,2,Domain, R,1,Range, Dune::FieldVector<Domain,1> > BasisTraits;
  public:
    typedef LocalFiniteElementTraits<
        C1LocalBasisVirtualInterface<BasisTraits>,
        LocalCoefficientsVirtualInterface,
        LocalInterpolationVirtualInterface< Domain, Range >
        > Traits;
    typedef typename Traits::LocalBasisType LocalBasis;
    typedef typename Traits::LocalCoefficientsType LocalCoefficients;
    typedef typename Traits::LocalInterpolationType LocalInterpolation;

    PQ22DLocalFiniteElement ( const GeometryType &gt )
      : gt_(gt),
        basis_(0), coefficients_(0), interpolation_(0)
    {
      if ( gt.isTriangle() )
        setup( Pk2DLocalFiniteElement<D,R,2>() );
      else if ( gt.isQuadrilateral() )
        setup( Q22DLocalFiniteElement<D,R>() );
    }
    PQ22DLocalFiniteElement ( const GeometryType &gt, const std::vector<unsigned int> vertexmap )
      : gt_(gt),
        basis_(0), coefficients_(0), interpolation_(0)
    {
      if ( gt.isTriangle() )
        setup( Pk2DLocalFiniteElement<D,R,2>(vertexmap) );
      else if ( gt.isQuadrilateral() )
        setup( Q22DLocalFiniteElement<D,R>() );
    }

    virtual ~PQ22DLocalFiniteElement ( )
    {
      delete coefficients_;
      delete interpolation_;
      delete basis_;
    }
    const LocalBasis& localBasis () const
    {
      return *basis_;
    }
    const LocalCoefficients& localCoefficients () const
    {
      return *coefficients_;
    }
    const LocalInterpolation& localInterpolation () const
    {
      return *interpolation_;
    }
    const GeometryType &type () const
    {
      return gt_;
    }
  private:
    template <class FE>
    void setup(const FE& fe)
    {
      basis_ = new C1LocalBasisVirtualImp< BasisTraits, typename FE::Traits::LocalBasisType >(fe);
      coefficients_ = new LocalCoefficientsVirtualImp< typename FE::Traits::LocalCoefficientsType >(fe);
      interpolation_ = new LocalInterpolationVirtualImp< Domain, Range, typename FE::Traits::LocalInterpolationType >(fe);
    }

    const GeometryType gt_;
    const LocalBasis *basis_;
    const LocalCoefficients *coefficients_;
    const LocalInterpolation *interpolation_;
  };

}

#endif
