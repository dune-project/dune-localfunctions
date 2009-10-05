// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH
#define DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH

#include <dune/grid/genericgeometry/referencemappings.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int dim, class QuadratureCreator >
    struct SubQuadratureCreator
    {
      typedef typename QuadratureCreator::Field Field;
      static const unsigned int dimension = dim;
      static const unsigned int codimension = dimension - QuadratureCreator::dimension;
      typedef GenericGeometry::Quadrature< dimension, Field > Quadrature;

      typedef std::pair< unsigned int, typename QuadratureCreator::Key > Key;

      const Quadrature &quadrature ( const unsigned int topologyyId, const Key &key )
      {
        typedef ReferenceMappings< Field, dimension > RefMappings;
        typedef typename RefMappings::Container RefMappingsContainer;
        typedef typename RefMappingsContainer::template Codim< codimension >::Mapping Mapping;

        const RefMappingsContainer &refMappings = RefMappings::container( topologyId );
        const unsigned int subTopologyId = refMappings.topology().topologyId( codimension, key.first );
        const Mapping &mapping = refMappings.template mapping< codimension >( key.first );

        Quadrature *subQuadrature = new Quadrature( topologyId );

        const typename QuadratureCreator::Quadrature *quadrature;
        IfTopology< QuadratureMaker, dimension >::apply( subTopologyId, key.second, quadrature );

        const unsigned int numQuadraturePoints = quadrature->size();
        for( unsigned int qp = 0; qp < numQuadraturePoints; ++qp )
          subQuadrature->insert( mapping.global( quadrature->point() ), quadrature->weight() );
        QuadratureCreator::release( *quadrature );

        return *subQuadrature;
      }

      template< class Topology >
      const Quadrature &quadrature ( const Key &key )
      {
        return quadrature( Topology::id, key );
      }

      static void release ( const Quadrature &quadrature )
      {
        delete &quadrature;
      }

    private:
      template< class Topology >
      struct QuadratureMaker
      {
        static void apply ( const typename QuadratureCreator::Key &key,
                            typename QuadratureCreator::Quadrature *&quadrature )
        {
          quadrature = QuadratureCreator::template quadrature< Topology >( key );
        }
      };
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_SUBQUADRATURE_HH
