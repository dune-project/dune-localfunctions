// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASISPROVIDER_HH
#define DUNE_BASISPROVIDER_HH

#include <map>
#include <dune/common/fvector.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // BasisProvider
  // -------------
  template< class BasisCreator >
  struct BasisProvider
  {
    typedef typename BasisCreator::StorageField StorageField;
    static const unsigned int dimension = BasisCreator::dimension;

    typedef typename BasisCreator::Key Key;
    typedef typename BasisCreator::Basis Basis;

    template< class Topology >
    static const Basis &basis ( const Key &key )
    {
      dune_static_assert( Topology::dimension == dimension, "basis provider called with a topology of wrong dimension!" );
      return instance().template getBasis< Topology >( key );
    }

    static const Basis &basis ( const unsigned int topologyId, const Key &key )
    {
      assert( topologyId < numTopologies );
      return instance().getBasis( topologyId, key );
    }

    static void release ( const Basis &basis )
    {}

  private:
    static const unsigned int numTopologies = (1 << dimension);

    typedef FieldVector< const Basis *, numTopologies > BasisArray;
    typedef std::map< Key, BasisArray > Storage;

    BasisProvider ()
    {}

    ~BasisProvider ()
    {
      const typename Storage::iterator end = basis_.end();
      for( typename Storage::iterator it = basis_.begin(); it != end; ++it )
      {
        for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        {
          const Basis *&basis = it->second[ topologyId ];
          if( basis != 0 )
            BasisCreator::release( *basis );
          basis = 0;
        }
      }
    }

    static BasisProvider &instance ()
    {
      static BasisProvider instance;
      return instance;
    }

    const Basis &getBasis ( const unsigned int topologyId, const Key &key )
    {
      typename Storage::iterator it = basis_.find( key );
      if( it == basis_.end() )
        it = basis_.insert( std::make_pair( key, BasisArray( 0 ) ) ).first;
      const Basis *&basis = it->second[ topologyId ];
      if( basis == 0 )
        GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, basis );
      return *basis;
    }

    template< class Topology >
    const Basis &getBasis ( const Key &key )
    {
      const unsigned int topologyId = Topology::id;
      typename Storage::iterator it = basis_.find( key );
      if( it == basis_.end() )
        it = basis_.insert( std::make_pair( key, BasisArray( 0 ) ) ).first;
      const Basis *&basis = it->second[ topologyId ];
      if( basis == 0 )
        basis = &BasisCreator::template basis< Topology >( key );
      return *basis;
    }

    template< class Topology >
    struct Maker
    {
      static void apply ( const Key &key, const Basis *&basis )
      {
        basis = &BasisCreator::template basis< Topology >( key );
      };
    };

    Storage basis_;
  };


  template< class BasisCreator, class CoefficientsCreator, class InterpolationCreator >
  struct FiniteElementProvider
  {
    typedef typename BasisCreator::StorageField StorageField;
    static const unsigned int dimension = BasisCreator::dimension;

    typedef typename BasisCreator::Key Key;
    typedef typename BasisCreator::Basis Basis;
    typedef typename CoefficientsCreator::LocalCoefficients Coefficients;
    typedef typename InterpolationCreator::LocalInterpolation Interpolation;

    struct FiniteElement
    {
      typedef typename BasisCreator::Basis Basis;
      typedef typename CoefficientsCreator::LocalCoefficients Coefficients;
      typedef typename InterpolationCreator::LocalInterpolation Interpolation;
      FiniteElement() : basis_(0), coeff_(0), interpol_(0) {}
      const Basis &basis() const
      {
        return *basis_;
      }
      const Coefficients &coefficients() const
      {
        return *coeff_;
      }
      const Interpolation &interpolation() const
      {
        return *interpol_;
      }
      void setBasis( const Basis *basis)
      {
        basis_ = basis;
      }
      void setCoefficients( const Coefficients *coeff)
      {
        coeff_ = coeff;
      }
      void setInterpolation( const Interpolation *interpol)
      {
        interpol_ = interpol;
      }
      bool empty() const
      {
        return basis_==0;
      }
      void release()
      {
        if (basis_)
          BasisCreator::release(*basis_);
        if (coeff_)
          CoefficientsCreator::release(*coeff_);
        if (interpol_)
          InterpolationCreator::release(*interpol_);
        basis_=0;
        coeff_=0;
        interpol_=0;
      }
    private:
      const Basis *basis_;
      const Coefficients *coeff_;
      const Interpolation *interpol_;
    };

    template< class Topology >
    static const FiniteElement &finiteElement ( const Key &key )
    {
      dune_static_assert( Topology::dimension == dimension, "basis provider called with a topology of wrong dimension!" );
      return instance().template getFiniteElement< Topology >( key );
    }
    static const FiniteElement &finiteElement ( const unsigned int topologyId, const Key &key )
    {
      assert( topologyId < numTopologies );
      return instance().getFiniteElement( topologyId, key );
    }

    static void release ( const FiniteElement & )
    {}

  private:
    static const unsigned int numTopologies = (1 << dimension);

    typedef FieldVector< FiniteElement, numTopologies > BasisArray;
    typedef std::map< Key, BasisArray > Storage;

    FiniteElementProvider ()
    {}

    ~FiniteElementProvider ()
    {
      const typename Storage::iterator end = fe_.end();
      for( typename Storage::iterator it = fe_.begin(); it != end; ++it )
      {
        for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        {
          FiniteElement &fe = it->second[ topologyId ];
          fe.release();
        }
      }
    }

    static FiniteElementProvider &instance ()
    {
      static FiniteElementProvider instance;
      return instance;
    }

    const FiniteElement &getFiniteElement ( const unsigned int topologyId, const Key &key )
    {
      typename Storage::iterator it = fe_.find( key );
      if( it == fe_.end() )
        it = fe_.insert( std::make_pair( key, BasisArray( ) ) ).first;
      FiniteElement &finiteElement = it->second[ topologyId ];
      if( finiteElement.empty() )
        GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, finiteElement );
      return finiteElement;
    }

    template< class Topology >
    const FiniteElement &getFiniteElement ( const Key &key )
    {
      const unsigned int topologyId = Topology::id;
      typename Storage::iterator it = fe_.find( key );
      if( it == fe_.end() )
        it = fe_.insert( std::make_pair( key, BasisArray( ) ) ).first;
      FiniteElement &finiteElement = it->second[ topologyId ];
      if( finiteElement.empty() )
        abort();
      return finiteElement;
    }

    template< class Topology >
    struct Maker
    {
      static void apply ( const Key &key, FiniteElement &finiteElement )
      {
        finiteElement.setBasis( &BasisCreator::template basis< Topology >( key ) );
        finiteElement.setCoefficients( &CoefficientsCreator::template localCoefficients< Topology >( key ) );
        finiteElement.setInterpolation( &InterpolationCreator::template localInterpolation< Topology >( key ) );
      };
    };

    Storage fe_;
  };

}

#endif // DUNE_BASISPROVIDER_HH
