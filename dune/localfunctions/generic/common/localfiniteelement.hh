// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/localfunctions/common/localfiniteelement.hh>

namespace Dune
{
  /**
   * @brief A LocalFiniteElement implementation bassed on three
   *        TopologyFactories providing the LocalBasis, LocalCoefficients,
   *        and LocalInterpolations. Note the key type for all three
   *        factories must coincide.
   **/
  template< class BasisF, class CoeffF, class InterpolF,
      unsigned int dimDomain, class D, class R >
  struct GenericLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits< typename BasisF::Object,
              typename CoeffF::Object,
              typename InterpolF::Object >,
          GenericLocalFiniteElement<BasisF, CoeffF, InterpolF, dimDomain,D,R> >
  {
    typedef GenericLocalFiniteElement<BasisF, CoeffF, InterpolF, dimDomain,D,R> This;
    typedef LocalFiniteElementInterface<
        LocalFiniteElementTraits<
            typename BasisF::Object,
            typename CoeffF::Object,
            typename InterpolF::Object >, This > Base;
    typedef typename Base::Traits Traits;

    typedef typename BasisF::Key Key;

    dune_static_assert( (Conversion<Key,typename CoeffF::Key>::sameType),
                        "incompatible keys between BasisCreator and CoefficientsCreator" );
    dune_static_assert( (Conversion<Key,typename InterpolF::Key>::sameType),
                        "incompatible keys between BasisCreator and InterpolationCreator" );

    /** \todo Please doc me !
     */
    GenericLocalFiniteElement ( unsigned int topologyId,
                                const Key &key )
      : topologyId_(topologyId),
        finiteElement_( )
    {
      GenericGeometry::IfTopology< FiniteElement::template Maker, dimDomain >::apply( topologyId, key, finiteElement_ );
    }
    ~GenericLocalFiniteElement()
    {
      finiteElement_.release();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return *(finiteElement_.basis_);
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return *(finiteElement_.coeff_);
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return *(finiteElement_.interpol_);
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      if ( GenericGeometry::hasGeometryType( topologyId_, dimDomain ) )
        return GenericGeometry::geometryType( topologyId_, dimDomain );
      return GeometryType();
    }

    /** \todo Please doc me !
     */
    unsigned int topologyId () const
    {
      return topologyId_;
    }
  private:
    struct FiniteElement
    {
      FiniteElement() : basis_(0), coeff_(0), interpol_(0) {}
      template <class Topology>
      void create( const Key &key )
      {
        release();
        basis_ = BasisF::template create<Topology>(key);
        coeff_ = CoeffF::template create<Topology>(key);
        interpol_ = InterpolF::template create<Topology>(key);
      }
      void release()
      {
        if (basis_)
          BasisF::release(basis_);
        if (coeff_)
          CoeffF::release(coeff_);
        if (interpol_)
          InterpolF::release(interpol_);
        basis_=0;
        coeff_=0;
        interpol_=0;
      }
      template< class Topology >
      struct Maker
      {
        static void apply ( const Key &key, FiniteElement &finiteElement )
        {
          finiteElement.template create<Topology>(key);
        };
      };
      typename Traits::LocalBasisType *basis_;
      typename Traits::LocalCoefficientsType *coeff_;
      typename Traits::LocalInterpolationType *interpol_;
    };
    unsigned int topologyId_;
    FiniteElement finiteElement_;
  };

}

#endif
