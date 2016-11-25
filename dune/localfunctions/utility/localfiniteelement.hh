// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/utility/l2interpolation.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>

namespace Dune
{
  /**
   * \brief A LocalFiniteElement implementation based on three
   *        TopologyFactories providing the LocalBasis, LocalCoefficients,
   *        and LocalInterpolations. Note the key type for all three
   *        factories must coincide.
   **/
  template< class BasisF, class CoeffF, class InterpolF>
  struct GenericLocalFiniteElement
  {
    typedef GenericLocalFiniteElement<BasisF, CoeffF, InterpolF> This;
    typedef LocalFiniteElementTraits< typename BasisF::Object,
        typename CoeffF::Object,
        typename InterpolF::Object > Traits;

    typedef typename BasisF::Key Key;
    static const unsigned int dimDomain = BasisF::dimension;

    typedef BasisF BasisFactory;
    typedef CoeffF CoefficientFactory;
    typedef InterpolF InterpolationFactory;

    static_assert(std::is_same<Key, typename CoeffF::Key>::value,
                  "incompatible keys between BasisCreator and CoefficientsCreator");
    static_assert(std::is_same<Key, typename InterpolF::Key>::value,
                  "incompatible keys between BasisCreator and InterpolationCreator" );

    /** \todo Please doc me */
    GenericLocalFiniteElement ( const GeometryType &gt, const Key &key )
      : topologyId_( gt.id() ),
        key_( key ),
        finiteElement_()
    {
      Impl::IfTopology< FiniteElement::template Maker, dimDomain >::apply( topologyId_, key_, finiteElement_ );
    }

    /** \todo Please doc me */
    GenericLocalFiniteElement ( const GenericLocalFiniteElement &other )
      : topologyId_( other.topologyId_ ),
        key_( other.key_ ),
        finiteElement_()
    {
      Impl::IfTopology< FiniteElement::template Maker, dimDomain >::apply( topologyId_, key_, finiteElement_ );
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

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return finiteElement_.basis_->size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return GeometryType(topologyId_,dimDomain);
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
    Key key_;
    FiniteElement finiteElement_;
  };

  /**
   * @brief Takes the basis and interpolation factory from a given
   *        LocalFiniteElement (derived from GenericLocalFiniteElement)
   *        and replaces the coefficients with dg local keys, i.e.,
   *        attaches all degrees of freedom to the codimension zero entity.
   **/
  template <class FE>
  struct DGLocalFiniteElement
    : public GenericLocalFiniteElement< typename FE::BasisFactory,
          DGLocalCoefficientsFactory< typename FE::BasisFactory >,
          typename FE::InterpolationFactory>
  {
    typedef GenericLocalFiniteElement< typename FE::BasisFactory,
        DGLocalCoefficientsFactory< typename FE::BasisFactory >,
        typename FE::InterpolationFactory> Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    DGLocalFiniteElement ( const GeometryType &gt, const typename Base::Key &key  )
      : Base( gt, key )
    {}
  };
  /**
   * @brief Takes the basis factory from a given
   *        LocalFiniteElement (derived from GenericLocalFiniteElement)
   *        and replaces the coefficients with dg local keys, i.e.,
   *        attaches all degrees of freedom to the codimension zero entity
   *        and uses a l2 interpolation.
   **/
  template <class FE>
  struct L2LocalFiniteElement
    : public GenericLocalFiniteElement< typename FE::BasisFactory,
          DGLocalCoefficientsFactory< typename FE::BasisFactory >,
          LocalL2InterpolationFactory< typename FE::BasisFactory, false > >
  {
    typedef GenericLocalFiniteElement< typename FE::BasisFactory,
        DGLocalCoefficientsFactory< typename FE::BasisFactory >,
        LocalL2InterpolationFactory< typename FE::BasisFactory, false > > Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    L2LocalFiniteElement ( const GeometryType &gt, const typename Base::Key &key  )
      : Base( gt, key )
    {}
  };
}

#endif
