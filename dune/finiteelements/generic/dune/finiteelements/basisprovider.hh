// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASISPROVIDER_HH
#define DUNE_BASISPROVIDER_HH
#include <fstream>
#include <map>
namespace Dune
{
#if 0
  template< class BasisCreator >
  struct BasisProvider
  {
    typedef typename BasisCreator :: StorageField StorageField;
    static const int dimension = BasisCreator :: dimension;
    typedef typename BasisCreator :: Basis Basis;

    static const Basis &basis(unsigned int id,unsigned int order)
    {
      return instance().getBasis(id,order);
    }
    template <class Topology>
    static const Basis &basis(unsigned int order)
    {
      return instance().template getBasis<Topology>(order);
    }
  private:
    enum { numTopologies = (1 << dimension) };
    BasisProvider()
    {}
    static BasisProvider &instance()
    {
      static BasisProvider instance;
      return instance;
    }
    const Basis &getBasis(unsigned int id,unsigned int order)
    {
      if (order>=basis_.size())
      {
        basis_.resize(order+1,FieldVector<Basis*,numTopologies>(0));
        GenericGeometry::IfTopology<BasisCreator::template Maker,dimension>::apply(id,order,basis_[order][id]);
      }
      else if (basis_[order][id] == 0)
        GenericGeometry::IfTopology<BasisCreator::template Maker,dimension>::apply(id,order,basis_[order][id]);
      return *(basis_[order][id]);
    }
    template <class Topology>
    const Basis &getBasis(unsigned int order)
    {
      const unsigned int id = Topology::id;
      if (order>=basis_.size())
      {
        basis_.resize(order+1,FieldVector<Basis*,numTopologies>(0));
        BasisCreator::template Maker<Topology>::apply(order,basis_[order][id]);
      }
      else if (basis_[order][id] == 0)
        BasisCreator::template Maker<Topology>::apply(order,basis_[order][id]);
      return *(basis_[order][id]);
    }
    std::vector<FieldVector<Basis*,numTopologies> > basis_;
  };
#endif
  template< class BasisCreator >
  struct BasisProvider
  {
    typedef typename BasisCreator :: StorageField StorageField;
    static const int dimension = BasisCreator :: dimension;
    typedef typename BasisCreator :: Basis Basis;
    typedef typename BasisCreator :: Key Key;
    static const Basis &basis(unsigned int id,Key key)
    {
      return instance().getBasis(id,key);
    }
    template <class Topology>
    static const Basis &basis(Key key)
    {
      return instance().template getBasis<Topology>(key);
    }

    static void release ( const Basis &basis )
    {}

  private:
    enum { numTopologies = (1 << dimension) };
    typedef std::map<Key,FieldVector<Basis*,numTopologies> > Storage;
    BasisProvider()
    {}
    static BasisProvider &instance()
    {
      static BasisProvider instance;
      return instance;
    }
    const Basis &getBasis(unsigned int id,Key key)
    {
      typename Storage::iterator b = basis_.find(key);
      if ( b == basis_.end() )
        b = basis_.insert(std::make_pair(key,FieldVector<Basis*,numTopologies>(0))).first;
      if (b->second[id] == 0)
        GenericGeometry::IfTopology<Maker,dimension>::apply(id,key,b->second[id]);
      return *(b->second[id]);
    }
    template <class Topology>
    const Basis &getBasis(Key key)
    {
      const unsigned int id = Topology::id;
      typename Storage::iterator b = basis_.find(key);
      if ( b == basis_.end() )
        b = basis_.insert(std::make_pair(key,FieldVector<Basis*,numTopologies>(0))).first;
      if (b->second[id] == 0)
        BasisCreator::template basis<Topology>(key,b->second[id]);
      return *(b->second[id]);
    }
    Storage basis_;
    template <class Topology>
    struct Maker {
      static void apply(Key key,Basis* &basis)
      {
        BasisCreator::template basis<Topology>(key,basis);
      };
    };
  };
}
#endif // DUNE_BASISPROVIDER_HH
