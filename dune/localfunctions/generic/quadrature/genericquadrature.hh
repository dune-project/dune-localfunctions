// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICQUADRATURE_HH
#define DUNE_GENERICQUADRATURE_HH

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/localfunctions/generic/quadrature/quadrature.hh>
#include <dune/localfunctions/generic/common/topologyfactory.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GenericQuadrature
    // -----------------

    /**
     * @brief extends a 1d quadrature to a generic reference elemenet
     *
     * \tparam Topology the topology of the reference element
     * \tparam OneDQuad a quadrature for [0,1]
     *
     * The 1d quadrature must be a std::vector-like container with a constructor
     * taking an order parameter
     **/
    template< class Topology, class OneDQuad >
    class GenericQuadrature;

    /** \brief Generic Quadrature for Point
    **/
    template< class OneDQuad >
    class GenericQuadrature< Point, OneDQuad >
      : public Quadrature< typename OneDQuad::Field, 0 >
    {
      typedef typename OneDQuad::Field F;
      typedef GenericQuadrature< Point, OneDQuad > This;
      typedef Quadrature< F, 0 > Base;

    public:
      typedef Point Topology;

      typedef F Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        Base::insert( Vector( Field( 0 ) ), 1 );
      }
    };


    /** \brief Generic Quadrature for Prisms
    **/
    template< class BaseTopology, class OneDQuad >
    class GenericQuadrature< Prism< BaseTopology >, OneDQuad >
      : public Quadrature< typename OneDQuad::Field, Prism< BaseTopology >::dimension >
    {
      typedef typename OneDQuad::Field F;
      typedef GenericQuadrature< Prism< BaseTopology >, OneDQuad > This;
      typedef Quadrature< F, Prism< BaseTopology >::dimension > Base;

    public:
      typedef Prism< BaseTopology > Topology;

      typedef F Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef OneDQuad OneDQuadrature;
      typedef GenericQuadrature< BaseTopology, OneDQuadrature > BaseQuadrature;

    public:
      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        OneDQuadrature onedQuad( order );
        BaseQuadrature baseQuad( order );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad[bqi].point( );
          const Field &baseWeight = baseQuad[bqi].weight( );

          Vector point;
          for( unsigned int i = 0; i < dimension-1; ++i )
            point[ i ] = basePoint[ i ];

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            point[ dimension-1 ] = onedQuad[oqi].point()[ 0 ];
            Base::insert( point, baseWeight * onedQuad[oqi].weight( ) );
          }
        }
      }
    };


    /** \brief Generic Quadrature for Pyramids
     *
     *  This quadrature for \f$B^\circ\f$ is generated from a quadrature for
     *  \f$B\f$ and a 1D quadrature by the so-called Duffy-Transformation
     *  \f$y(x,z) = ((1-z)x,y)^T\f$. Hence, we have
     *  \f[
     *  \int_{B^\circ} f( y )\,\mathrm{d}y
     *  = \int_0^1 \int_B f( (1-z)x, z )\,\mathrm{d}x\,(1-z)^{\dim B}\,\mathrm{d}z.
     *  \f]
     *  Therefore, the 1D quadrature must be at least \f$\dim B\f$ orders higher
     *  than the quadrature for \f$B\f$.
     *
     *  Question: If the polynomials are created via Duffy Transformation, do we
     *            really need a higher quadrature order?
     */
    template< class BaseTopology, class OneDQuad >
    class GenericQuadrature< Pyramid< BaseTopology >, OneDQuad >
      : public Quadrature< typename OneDQuad::Field, Pyramid< BaseTopology >::dimension >
    {
      typedef typename OneDQuad::Field F;
      typedef GenericQuadrature< Pyramid< BaseTopology >, OneDQuad > This;
      typedef Quadrature< F, Pyramid< BaseTopology >::dimension > Base;

    public:
      typedef Pyramid< BaseTopology > Topology;

      typedef F Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef OneDQuad OneDQuadrature;
      typedef GenericQuadrature< BaseTopology, OneDQuadrature > BaseQuadrature;

    public:
      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        OneDQuadrature onedQuad( order + dimension-1 );
        BaseQuadrature baseQuad( order );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad[bqi].point( );
          const Field &baseWeight = baseQuad[bqi].weight( );

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            Vector point;
            point[ dimension-1 ] = onedQuad[oqi].point( )[ 0 ];
            const Field scale = Field( 1 ) - point[ dimension-1 ];
            for( unsigned int i = 0; i < dimension-1; ++i )
              point[ i ] = scale * basePoint[ i ];

            Field weight = baseWeight * onedQuad[oqi].weight( );
            for ( unsigned int p = 0; p<dimension-1; ++p)
              weight *= scale;                    // pow( scale, dimension-1 );
            Base::insert( point, weight );
          }
        }
      }
    };


    /**
     * @brief Factory for the generic quadratures
     *
     * This is a Dune::GenericGeometry::TopologyFactory creating
     * GenericQuadrature from a given 1d quadrature
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     * \tparam OneDQuad the underlying 1d quadrature
     *
     * Note: the computation of the quadrature points and weights are
     * carried out in the field type of the 1d quadrature which can differ from F.
     **/
    template< int dim, class F, class OneDQuad >
    struct GenericQuadratureFactory;

    template< int dim, class F, class OneDQuad >
    struct GenericQuadratureFactoryTraits
    {
      static const unsigned int dimension = dim;
      typedef unsigned int Key;
      typedef const Quadrature<F,dim> Object;
      typedef GenericQuadratureFactory<dim,F,OneDQuad> Factory;
    };

    template< int dim, class F, class OneDQuad >
    struct GenericQuadratureFactory :
      public TopologyFactory< GenericQuadratureFactoryTraits<dim,F,OneDQuad> >
    {
      static const unsigned int dimension = dim;
      typedef F Field;
      typedef GenericQuadratureFactoryTraits<dim,F,OneDQuad> Traits;

      typedef typename Traits::Key Key;
      typedef typename Traits::Object Object;

      template< class Topology >
      static Object* createObject ( const Key &order )
      {
        return new Object( GenericQuadrature< Topology, OneDQuad >( order ) );
      }
    };


    // GenericQuadratureProvider
    // ---------------------

    /**
     * @brief Singleton factory for the generic quadratures
     *
     * Wrapper for the Dune::GenericGeometry::GenericQuadratureFactory providing
     * singleton storage.
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     * \tparam OneDQuad the underlying 1d quadrature
     **/
    template< int dim, class F, class OneDQuad >
    struct GenericQuadratureProvider
      : public TopologySingletonFactory< GenericQuadratureFactory< dim, F, OneDQuad > >
    {};


  }

}

#endif
