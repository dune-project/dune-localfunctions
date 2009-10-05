// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICQUADRATURE_HH
#define DUNE_GENERICQUADRATURE_HH

#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>
#include <dune/finiteelements/quadrature/gaussquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // DefaultQuadratureTraits
    // -----------------------

    struct DefaultQuadratureTraits
    {
      typedef double Field;

      typedef GaussQuadrature< Field > OneDQuadrature;
    };



    // GenericQuadrature
    // -----------------

    template< class Topology, class Traits = DefaultQuadratureTraits >
    class GenericQuadrature;


    template< class Traits >
    class GenericQuadrature< Point, Traits >
      : public Quadrature< 0, typename Traits::Field >
    {
      typedef GenericQuadrature< Point, Traits > This;
      typedef Quadrature< 0, typename Traits::Field > Base;

    public:
      typedef Point Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        Base::insert( Vector( Field( 0 ) ), 1 );
      }
    };


    template< class BaseTopology, class Traits >
    class GenericQuadrature< Prism< BaseTopology >, Traits >
      : public Quadrature< Prism< BaseTopology >::dimension, typename Traits::Field >
    {
      typedef GenericQuadrature< Prism< BaseTopology >, Traits > This;
      typedef Quadrature< Prism< BaseTopology >::dimension, typename Traits::Field > Base;

    public:
      typedef Prism< BaseTopology > Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef typename Traits::OneDQuadrature OneDQuadrature;
      typedef GenericQuadrature< BaseTopology, Traits > BaseQuadrature;

    public:
      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        OneDQuadrature onedQuad( order );
        BaseQuadrature baseQuad( order );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad.point( bqi );
          const Field &baseWeight = baseQuad.weight( bqi );

          Vector point;
          for( unsigned int i = 0; i < dimension-1; ++i )
            point[ i ] = basePoint[ i ];

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            const typename OneDQuadrature::Vector onedPoint = onedQuad.point( oqi );
            const Field &onedWeight = onedQuad.weight( oqi );

            point[ dimension-1 ] = onedPoint[ 0 ];
            Base::insert( point, baseWeight * onedWeight );
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
    template< class BaseTopology, class Traits >
    class GenericQuadrature< Pyramid< BaseTopology >, Traits >
      : public Quadrature< Pyramid< BaseTopology >::dimension, typename Traits::Field >
    {
      typedef GenericQuadrature< Pyramid< BaseTopology >, Traits > This;
      typedef Quadrature< Pyramid< BaseTopology >::dimension, typename Traits::Field > Base;

    public:
      typedef Pyramid< BaseTopology > Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef typename Traits::OneDQuadrature OneDQuadrature;
      typedef GenericQuadrature< BaseTopology, Traits > BaseQuadrature;

    public:
      explicit GenericQuadrature ( unsigned int order )
        : Base( Topology::id )
      {
        OneDQuadrature onedQuad( order + dimension-1 );
        BaseQuadrature baseQuad( order );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad.point( bqi );
          const Field &baseWeight = baseQuad.weight( bqi );

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            const typename OneDQuadrature::Vector onedPoint = onedQuad.point( oqi );
            const Field &onedWeight = onedQuad.weight( oqi );

            Vector point;
            point[ dimension-1 ] = onedPoint[ 0 ];
            const Field scale = Field( 1 ) - onedPoint[ 0 ];
            for( unsigned int i = 0; i < dimension-1; ++i )
              point[ i ] = scale * basePoint[ i ];

            const Field weight = baseWeight * onedWeight * pow( scale, dimension-1 );
            Base::insert( point, weight );
          }
        }
      }
    };

  }

}

#endif
