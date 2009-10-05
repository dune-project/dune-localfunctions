// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH
#define DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH

#include <dune/grid/genericgeometry/referencemappings.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // SubQuadratureCreator
    // --------------------

    template< unsigned int dim, class QuadratureCreator >
    struct SubQuadratureCreator
    {
      typedef typename QuadratureCreator::Field Field;
      static const unsigned int dimension = dim;
      static const unsigned int codimension = dimension - QuadratureCreator::dimension;
      typedef GenericGeometry::Quadrature< dimension, Field > Quadrature;
      typedef typename QuadratureCreator::Quadrature SubQuadrature;

      typedef std::pair< unsigned int, typename QuadratureCreator::Key > Key;

      static const Quadrature &
      quadrature ( const unsigned int topologyId, const Key &key )
      {
        typedef ReferenceMappings< Field, dimension > RefMappings;
        typedef typename RefMappings::Container RefMappingsContainer;
        typedef typename RefMappingsContainer::template Codim< codimension >::Mapping Mapping;

        const RefMappingsContainer &refMappings = RefMappings::container( topologyId );
        const unsigned int subTopologyId = refMappings.topology().topologyId( codimension, key.first );
        const Mapping &mapping = refMappings.template mapping< codimension >( key.first );

        Quadrature *subQuadrature = new Quadrature( topologyId );

        const typename QuadratureCreator::Quadrature *quadrature;
        IfTopology< QuadratureMaker, dimension-codimension >::apply( subTopologyId, key.second, quadrature );

        const unsigned int numQuadraturePoints = quadrature->size();
        for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
          subQuadrature->insert( mapping.global( quadrature->point( qp ) ), quadrature->weight( qp ) );
        QuadratureCreator::release( *quadrature );

        return *subQuadrature;
      }

      static const SubQuadrature &
      subQuadrature ( const unsigned int topologyId, const Key &key )
      {
        const unsigned int subTopologyId =
          ReferenceTopologies<dimension>::get(topologyId).topologyId( codimension,key.first );
        const SubQuadrature *subQuadrature;
        IfTopology< QuadratureMaker, dimension-codimension >::apply( subTopologyId, key.second, subQuadrature );
        return *subQuadrature;
      }

      template< class Topology >
      static const Quadrature &quadrature ( const Key &key )
      {
        return quadrature( Topology::id, key );
      }

      static void release ( const Quadrature &quadrature )
      {
        delete &quadrature;
      }

      void release( const SubQuadrature &quad )
      {
        QuadratureCreator::release( quad );
      }

    private:
      template< class Topology >
      struct QuadratureMaker
      {
        static void apply ( const typename QuadratureCreator::Key &key,
                            const typename QuadratureCreator::Quadrature *&quadrature )
        {
          quadrature = &(QuadratureCreator::template quadrature< Topology >( key ));
        }
      };
    };



    // SubQuadratureProvider
    // ---------------------

    template< unsigned int dim, class QuadratureCreator >
    struct SubQuadratureProvider
      : public QuadratureProvider< SubQuadratureCreator< dim, QuadratureCreator > >
    {
      typedef SubQuadratureCreator< dim, QuadratureCreator> Creator;
      typedef typename QuadratureCreator::Quadrature SubQuadrature;
      static const SubQuadrature &
      subQuadrature(unsigned int topologyId, const typename Creator::Key &key)
      {
        return Creator::subQuadrature(topologyId,key);
      }
      template <class Topology>
      static const SubQuadrature &
      subQuadrature(const typename Creator::Key &key)
      {
        return subQuadrature(Topology::id,key);
      }
      void release( const SubQuadrature &quad )
      {
        Creator::release( quad );
      }
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH
