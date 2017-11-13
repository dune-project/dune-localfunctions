// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQ22DLOCALFINITEELEMENT_HH
#define DUNE_PQ22DLOCALFINITEELEMENT_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>
#include "qk.hh"
#include "pk2d.hh"

namespace Dune
{
  template<class D, class R>
  class PQ22DLocalFiniteElement
  {
    typedef Dune::FieldVector<D,2> Domain;
    typedef Dune::FieldVector<R,1> Range;
    typedef LocalBasisTraits<D,2,Domain, R,1,Range, Dune::FieldMatrix<R,1,2>, 0 > BasisTraits;

    typedef typename Dune::LocalFiniteElementVirtualInterface<BasisTraits> LocalFiniteElementBase;
  public:
    typedef LocalFiniteElementTraits<
        LocalBasisVirtualInterface<BasisTraits>,
        LocalCoefficientsVirtualInterface,
        LocalInterpolationVirtualInterface< Domain, Range >
        > Traits;
    typedef typename Traits::LocalBasisType LocalBasis;
    typedef typename Traits::LocalCoefficientsType LocalCoefficients;
    typedef typename Traits::LocalInterpolationType LocalInterpolation;

    PQ22DLocalFiniteElement ( const GeometryType &gt )
      : gt_(gt)
    {
      if ( gt.isTriangle() )
        setup( Pk2DLocalFiniteElement<D,R,2>() );
      else if ( gt.isQuadrilateral() )
        setup( QkLocalFiniteElement<D,R,2,2>() );
    }

    PQ22DLocalFiniteElement ( const GeometryType &gt, const std::vector<unsigned int> vertexmap )
      : gt_(gt)
    {
      if ( gt.isTriangle() )
        setup( Pk2DLocalFiniteElement<D,R,2>(vertexmap) );
      else if ( gt.isQuadrilateral() )
        setup( QkLocalFiniteElement<D,R,2,2>() );
    }

    PQ22DLocalFiniteElement ( const PQ22DLocalFiniteElement<D, R>& other )
      : gt_(other.gt_)
    {
      fe_ = other.fe_->clone();
    }

    ~PQ22DLocalFiniteElement ( )
    {
      delete fe_;
    }

    const LocalBasis& localBasis () const
    {
      return fe_->localBasis();
    }

    const LocalCoefficients& localCoefficients () const
    {
      return fe_->localCoefficients();
    }

    const LocalInterpolation& localInterpolation () const
    {
      return fe_->localInterpolation();
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return fe_->localBasis().size();
    }

    const GeometryType &type () const
    {
      return gt_;
    }

  private:

    template <class FE>
    void setup(const FE& fe)
    {
      fe_ = new LocalFiniteElementVirtualImp<FE>(fe);
    }

    const GeometryType gt_;
    const LocalFiniteElementBase *fe_;
  };

}

#endif
