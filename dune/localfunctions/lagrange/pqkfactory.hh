// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQK_FACTORY_HH
#define DUNE_PQK_FACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>

namespace Dune
{

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Empty default implementation
   */
  template<class D, class R, int d, int k>
  struct DimSpecificPQkLocalFiniteElementFactory
  {
    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,d>::Traits::LocalBasisType::Traits,0>::Traits T;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      return 0;
    }
  };

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Specialization for dim=3
   */
  template<class D, class R, int k>
  struct DimSpecificPQkLocalFiniteElementFactory<D,R,3,k>
  {
    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,3>::Traits::LocalBasisType::Traits,0>::Traits T;
    typedef PrismP1LocalFiniteElement<D,R> PrismP1;
    typedef PrismP2LocalFiniteElement<D,R> PrismP2;
    typedef PyramidP1LocalFiniteElement<D,R> PyramidP1;
    typedef PyramidP2LocalFiniteElement<D,R> PyramidP2;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      if ((gt.isPrism())and (k==1))
        return new LocalFiniteElementVirtualImp<PrismP1>(PrismP1());
      if ((gt.isPrism())and (k==2))
        return new LocalFiniteElementVirtualImp<PrismP2>(PrismP2());
      if ((gt.isPyramid())and (k==1))
        return new LocalFiniteElementVirtualImp<PyramidP1>(PyramidP1());
      if ((gt.isPyramid())and (k==2))
        return new LocalFiniteElementVirtualImp<PyramidP2>(PyramidP2());
      return 0;
    }
  };


  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int k>
  struct PQkLocalFiniteElementFactory
  {
    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,0>::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    typedef P0LocalFiniteElement<D,R,dim> P0;
    typedef PkLocalFiniteElement<D,R,dim,k> Pk;
    typedef QkLocalFiniteElement<D,R,dim,k> Qk;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      if (k==0)
        return new LocalFiniteElementVirtualImp<P0>(P0(gt));

      if (gt.isSimplex())
        return new LocalFiniteElementVirtualImp<Pk>(Pk());

      if (gt.isCube())
        return new LocalFiniteElementVirtualImp<Qk>(Qk());

      return DimSpecificPQkLocalFiniteElementFactory<D,R,dim,k>::create(gt);
    }
  };



  /** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
   *
   * An interface for dealing with different vertex orders is currently missing.
   * So you can in general only use this for order=1,2 or with global DG spaces
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for shape function values
   * \tparam dim Element dimension
   * \tparam k Element order
   */
  template<class D, class R, int dim, int k>
  class PQkLocalFiniteElementCache
  {
  protected:
    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,0>::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PQkLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    PQkLocalFiniteElementCache(const PQkLocalFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~PQkLocalFiniteElementCache()
    {
      typename FEMap::iterator it = cache_.begin();
      typename FEMap::iterator end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt) const
    {
      typename FEMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        FiniteElementType* fe = PQkLocalFiniteElementFactory<D,R,dim,k>::create(gt);
        if (fe==0)
          DUNE_THROW(Dune::NotImplemented,"No Pk/Qk like local finite element available for geometry type " << gt << " and order " << k);

        cache_[gt] = fe;
        return *fe;
      }
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

}

#endif
