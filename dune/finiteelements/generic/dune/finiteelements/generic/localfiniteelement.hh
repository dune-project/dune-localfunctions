// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/generic/basisprovider.hh>

namespace Dune
{
  // forward declaration
  template< unsigned int dimDomain, class D, class R,
      class Basis >
  class GenericLocalBasis;

  template< class Interpolation >
  class GenericLocalInterpolation;

  template< class BasisC, class CoeffC, class InterpolC,
      unsigned int dimDomain, class D, class R >
  struct GenericLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename BasisC::Basis>,
              typename CoeffC::LocalCoefficients,
              GenericLocalInterpolation<typename InterpolC::LocalInterpolation> >,
          GenericLocalFiniteElement<BasisC, CoeffC, InterpolC, dimDomain,D,R> >
  {
    typedef FiniteElementProvider<BasisC,CoeffC,InterpolC> FECreator;
    typedef typename FECreator::FiniteElement FiniteElement;

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename FECreator::Basis>,
        typename FECreator::LocalCoefficients,
        GenericLocalInterpolation<typename FECreator::LocalInterpolation> > Traits;

    /** \todo Please doc me !
     */
    GenericLocalFiniteElement ( unsigned int topologyId,
                                unsigned int order )
      : topologyId_(topologyId),
        order_(order),
        finiteElement_( FECreator::finiteElement(topologyId,order) ),
        localBasis_(finiteElement_.basis()),
        localInterpolation_(finiteElement_.interpolation())
    {}
    ~GenericLocalFiniteElement()
    {
      FECreator::release( finiteElement_ );
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return localBasis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return finiteElement_.coefficients();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return localInterpolation_;
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
    unsigned int topologyId_;
    unsigned int order_;
    const FiniteElement &finiteElement_;
    typename Traits::LocalBasisType localBasis_;
    typename Traits::LocalInterpolationType localInterpolation_;
  };

  template< unsigned int dimDomain, class D, class R,
      class Basis >
  class GenericLocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,dimDomain,FieldVector<D,dimDomain>,R,1,FieldVector<R,1>,
            FieldVector<FieldVector<R,dimDomain>,1> >,
        GenericLocalBasis<dimDomain,D,R,Basis >
        >
  {
  public:
#if 0
    /** \brief Export the number of degrees of freedom */
    enum {N = (k+1)*(k+2)/2};
    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};
#endif
    typedef C1LocalBasisTraits<D,dimDomain,FieldVector<D,dimDomain>,R,1,FieldVector<R,1>,
        FieldVector<FieldVector<R,dimDomain>,1> > Traits;

    //! \brief Standard constructor
    GenericLocalBasis (const Basis &basis)
      : basis_(basis)
    {}

    //! \brief number of shape functions
    unsigned int size () const
    {
      return basis_.size();
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      basis_.evaluate(x,out);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,         // position
                      std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(size());
      basis_.jacobian(x,out);
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return basis_.order();
    }

  private:
    const Basis &basis_;
  };

  template< class Interpolation >
  struct GenericLocalInterpolation
    : public LocalInterpolationInterface< Interpolation >
  {
    GenericLocalInterpolation( const Interpolation &interpol)
      : interpol_(interpol)
    {}
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      interpol_.interpolate(f,out);
    }
  private:
    const Interpolation &interpol_;
  };
}

#endif
