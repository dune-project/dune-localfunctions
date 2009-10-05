// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURE_HH
#define DUNE_QUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // QuadraturePoint
    // ---------------

    template< class F, unsigned int dim >
    class QuadraturePoint
    {
      typedef QuadraturePoint< F, dim > This;

    public:
      static const unsigned int dimension = dim;

      typedef F Field;

      typedef FieldVector< Field, dimension > Vector;

    private:
      Vector position_;
      Field weight_;

    public:
      QuadraturePoint ( const Vector &position, const Field &weight )
        : position_( position ),
          weight_( weight )
      {}

      const Vector &position () const
      {
        return position_;
      }

      const Field &weight () const
      {
        return weight_;
      }
    };



    // QuadratureRule
    // --------------

    template< unsigned int dim, class F >
    class QuadratureRule
    {
      typedef QuadratureRule< dim, F > This;

    public:
      typedef F Field;
      static const unsigned int dimension = dim;

      typedef GenericGeometry::QuadraturePoint< Field, dimension > QuadraturePoint;
      typedef typename QuadraturePoint::Vector Vector;

    private:
      typedef std::vector< QuadraturePoint > QuadraturePointContainer;

    public:
      typedef typename QuadraturePointContainer::const_iterator iterator;

    private:
      std::vector< QuadraturePoint > points_;
      GeometryType type_;

    protected:
      explicit QuadratureRule ( const GeometryType &type )
        : type_( type )
      {}

    public:
      iterator begin () const
      {
        return points_.begin();
      }

      iterator end () const
      {
        return points_.end();
      }

      GeometryType type () const
      {
        return type_;
      }

      size_t size () const
      {
        return points_.size();
      }

    protected:
      void insert ( const QuadraturePoint &point )
      {
        points_.push_back( point );
      }

      void insert ( const Vector &position, const Field &weight )
      {
        insert( QuadraturePoint( position, weight ) );
      }
    };

  }

}

#endif
