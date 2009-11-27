// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURE_HH
#define DUNE_QUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/finiteelements/generic/math/field.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // QuadraturePoint
    // ---------------

    /**
     * @brief Base class for a single quadrature point (point and weight)
     *        with template field type and dimension
     **/
    template< class F, unsigned int dim >
    class QuadraturePoint
    {
      typedef QuadraturePoint< F, dim > This;

    public:
      static const unsigned int dimension = dim;

      typedef F Field;

      typedef FieldVector< Field, dimension > Vector;

    private:
      Vector point_;
      Field weight_;

    public:
      QuadraturePoint ( const Vector &point, const Field &weight )
        : point_( point ),
          weight_( weight )
      {}

      template <class FF>
      QuadraturePoint ( const QuadraturePoint<FF,dim>& other )
        : point_( field_cast<F>(other.point() ) ),
          weight_( field_cast<F>(other.weight() ) )
      {}

      const Vector &point () const
      {
        return point_;
      }

      const Field &weight () const
      {
        return weight_;
      }

    };


    // Quadrature
    // ----------

    /**
     * @brief Base class for the generic quadrature implementations.
     *        Field type and dimension are
     *        template argument construction is through a topology id.
     **/
    template< class F, unsigned int dim >
    class Quadrature
    {
      typedef Quadrature< F, dim > This;

    public:
      typedef F Field;
      static const unsigned int dimension = dim;

      typedef GenericGeometry::QuadraturePoint< Field, dimension > QuadraturePoint;
      typedef typename QuadraturePoint::Vector Vector;

      typedef typename std::vector< QuadraturePoint >::const_iterator Iterator;

    public:
      //! Constructor taking topology id
      explicit Quadrature ( const unsigned int topologyId )
        : topologyId_( topologyId )
      {}

      //! Copy constructor
      template< class Q >
      Quadrature ( const Q &q )
        : topologyId_( q.topologyId() )
      {
        points_.reserve( q.size() );
        const typename Q::Iterator end = q.end();
        for( typename Q::Iterator it = q.begin(); it != end; ++it )
          points_.push_back( *it );
      }

      //! Access a quadrature point
      const QuadraturePoint &operator[] ( const unsigned int i ) const
      {
        return points_[ i ];
      }

      //! start iterator over the quadrature points
      Iterator begin () const
      {
        return points_.begin();
      }

      //! end iterator over the quadrature points
      Iterator end () const
      {
        return points_.end();
      }

      //! access the coordinates of a quadrature point
      const Vector &point ( const unsigned int i ) const
      {
        return (*this)[ i ].point();
      }

      //! access the weight of a quadrature point
      const Field &weight ( const unsigned int i ) const
      {
        return (*this)[ i ].weight();
      }

      //! topology id of the quadrature
      unsigned int topologyId () const
      {
        return topologyId_;
      }

      //! number of quadrature points
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

      template< unsigned int d, class QC >
      friend struct SubQuadratureCreator;

    private:
      std::vector< QuadraturePoint > points_;
      unsigned int topologyId_;
    };

  }
}

#endif
