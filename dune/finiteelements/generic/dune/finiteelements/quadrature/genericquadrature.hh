// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICQUADRATURE_HH
#define DUNE_GENERICQUADRATURE_HH

#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/grid/quadrature/quadrature.hh>
#include <dune/grid/quadrature/gaussquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // DefaultQuadratureTraits
    // -----------------------

    struct DefaultQuadratureTraits
    {
      typedef double Field;

      static const GeometryType::BasicType linetype = GeometryType::simplex;

      typedef GaussQuadratureRule< Field > OneDQuadrature;
    };



    // GenericQuadratureRule
    // ---------------------

    template< class Topology, class Traits = DefaultQuadratureTraits >
    class GenericQuadratureRule;


    template< class Traits >
    class GenericQuadratureRule< Point, Traits >
      : public QuadratureRule< 0, typename Traits::Field >
    {
      typedef GenericQuadratureRule< Point, Traits > This;
      typedef QuadratureRule< 0, typename Traits::Field > Base;

    public:
      typedef Point Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

      explicit GenericQuadratureRule ( unsigned int order )
        : Base( DuneGeometryType< Topology, Traits::linetype >::type() )
      {
        Base::insert( Vector( Field( 0 ) ), 1 );
      }
    };


    template< class BaseTopology, class Traits >
    class GenericQuadratureRule< Prism< BaseTopology >, Traits >
      : public QuadratureRule< Prism< BaseTopology >::dimension, typename Traits::Field >
    {
      typedef GenericQuadratureRule< Prism< BaseTopology >, Traits > This;
      typedef QuadratureRule< Prism< BaseTopology >::dimension, typename Traits::Field > Base;

    public:
      typedef Prism< BaseTopology > Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef typename Traits::OneDQuadrature OneDQuadrature;
      typedef GenericQuadratureRule< BaseTopology, Traits > BaseQuadrature;

    public:
      explicit GenericQuadratureRule ( unsigned int order )
        : Base( DuneGeometryType< Topology, Traits::linetype >::type() )
      {
        typedef typename OneDQuadrature::iterator OneDQuadIterator;
        typedef typename BaseQuadrature::iterator BaseQuadIterator;

        OneDQuadrature onedQuad( order );
        BaseQuadrature baseQuad( order );

        const BaseQuadIterator bend = baseQuad.end();
        for( BaseQuadIterator bit = baseQuad.begin(); bit != bend; ++bit )
        {
          const typename BaseQuadrature::Vector &basePoint = bit->position();
          const Field &baseWeight = bit->weight();

          Vector point;
          for( unsigned int i = 0; i < dimension-1; ++i )
            point[ i ] = basePoint[ i ];

          const OneDQuadIterator oend = onedQuad.end();
          for( OneDQuadIterator oit = onedQuad.begin(); oit != oend; ++oit )
          {
            point[ dimension-1 ] = oit->position()[ 0 ];
            Base::insert( point, baseWeight * oit->weight() );
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
    class GenericQuadratureRule< Pyramid< BaseTopology >, Traits >
      : public QuadratureRule< Pyramid< BaseTopology >::dimension, typename Traits::Field >
    {
      typedef GenericQuadratureRule< Pyramid< BaseTopology >, Traits > This;
      typedef QuadratureRule< Pyramid< BaseTopology >::dimension, typename Traits::Field > Base;

    public:
      typedef Pyramid< BaseTopology > Topology;

      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Vector Vector;

    private:
      typedef typename Traits::OneDQuadrature OneDQuadrature;
      typedef GenericQuadratureRule< BaseTopology, Traits > BaseQuadrature;

    public:
      explicit GenericQuadratureRule ( unsigned int order )
        : Base( DuneGeometryType< Topology, Traits::linetype >::type() )
      {
        typedef typename OneDQuadrature::iterator OneDQuadIterator;
        typedef typename BaseQuadrature::iterator BaseQuadIterator;

        OneDQuadrature onedQuad( order + dimension-1 );
        BaseQuadrature baseQuad( order );

        const BaseQuadIterator bend = baseQuad.end();
        for( BaseQuadIterator bit = baseQuad.begin(); bit != bend; ++bit )
        {
          const typename BaseQuadrature::Vector &basePoint = bit->position();
          const Field &baseWeight = bit->weight();

          const OneDQuadIterator oend = onedQuad.end();
          for( OneDQuadIterator oit = onedQuad.begin(); oit != oend; ++oit )
          {
            Vector point;
            point[ dimension-1 ] = oit->position()[ 0 ];
            const Field scale = Field( 1 ) - point[ dimension-1 ];
            for( unsigned int i = 0; i < dimension-1; ++i )
              point[ i ] = scale * basePoint[ i ];

            const Field weight = baseWeight * oit->weight() * pow( scale, dimension-1 );
            Base::insert( point, weight );
          }
        }
      }
    };

  }

}

#endif
