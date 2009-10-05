// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/genericgeometry/referencemappings.hh>
namespace Dune {
  template <class Field>
  struct LabattoPoints1D
  {
    LabattoPoints1D(unsigned int order)
    {
      points_.resize(order-1);
      for (unsigned int i=1; i<order; ++i)
        points_[i-1] = Field(i)/Field(order);
      if (order==3) {
        points_[0] = 0.1;
        points_[1] = 0.9;
      }
      if (order==4) {
        points_[0] = 0.15;
        points_[1] = 0.5;
        points_[2] = 0.85;
      }
      if (order==5) {
        points_[0] = 0.1;
        points_[1] = 0.3;
        points_[2] = 0.7;
        points_[3] = 0.9;
      }
      if (order==6) {
        points_[0] = 0.09;
        points_[1] = 0.28;
        points_[2] = 0.5;
        points_[3] = 0.72;
        points_[4] = 0.91;
      }
    }
    const unsigned int size()
    {
      return points_.size();
    }
    const unsigned int order()
    {
      return points_.size()+1;
    }
    const Field &point(int i)
    {
      return points_[i];
    }
    std::vector<Field> points_;
  };

  template <class Field,class Topology>
  struct LabattoInnerPoints;

  template <class Field>
  struct LabattoInnerPoints<Field,GenericGeometry::Point>
  {
    static const unsigned int dimension = 0;
    static unsigned int size(const unsigned int order)
    {
      return 1;
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int i = order+1-iup;
      for ( unsigned int d=0; d<dimStart; ++d )
      {
        points[startPos].point_[d] -= points1D[i-2];
      }

      return startPos+1;
    }
    template <int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      points->point_[0] = 0;
    }
  };
  template <class Field,class Base>
  struct LabattoInnerPoints<Field,GenericGeometry::Pyramid<Base> >
  {
    typedef LabattoInnerPoints<Field,Base> LabattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      if (order<=dimension) return 0;
      unsigned int size=0;
      for (unsigned int o=1; o<=order-dimension; ++o)
        size += LabattoBase::size(o+1);
      return size;
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int endPos;
      for (unsigned int i=2; i<=order-iup; ++i)
      {
        endPos = LabattoBase::template setupSimplex<dim>(iup+i-1,dimStart,startPos,points1D,points);
        for (unsigned int j=startPos; j<endPos; ++j)
        {
          for (unsigned int d=0; d<dimStart; ++d)
          {
            if ( d==dimStart-dimension )
              points[j].point_[d] += dimStart*points1D[i-2];
            else
              points[j].point_[d] -= points1D[i-2];
          }
        }
        startPos = endPos;
      }
      return endPos;
    }
    template <unsigned int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int startPos=0,endPos=0;
      for (unsigned int i=2; i<=order; ++i)
      {
        endPos = LabattoBase::template setupSimplex<dim>(i-1,dimension,startPos,points1D,points);

        for (unsigned int j=startPos; j<endPos; ++j)
        {
          for (unsigned int d=0; d<dimension; ++d)
          {
            if ( d==0 )
              points[j].point_[d] += 1.+dimension*points1D[i-2];
            else
              points[j].point_[d] += 1.-points1D[i-2];
            points[j].point_[d] /= (dimension+1);
          }
        }
        startPos = endPos;
      }
    }
  };
  template <class Field,class Base>
  struct LabattoInnerPoints<Field,GenericGeometry::Prism<Base> >
  {
    typedef LabattoInnerPoints<Field,Base> LabattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      return LabattoBase::size(order)*(order-1);
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      abort();
    }
    template <unsigned int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      assert(dim>=dimension);
      assert(points1D.size()==order-1);
      LabattoBase::template setup<dim>(points1D,points);
      const unsigned int baseSize = LabattoBase::size(order);
      for (unsigned int q=0; q<order-1; q++)
      {
        for (unsigned int i=0; i<baseSize; ++i)
        {
          const unsigned int pos = q*baseSize+i;
          for (unsigned int d=0; d<dimension-1; ++d)
            points[pos].point_[d] = points[i].point_[d];
          points[pos].point_[dimension-1]=points1D[q];
        }
      }
    }
  };

  template <class Field,int dim>
  struct LabattoPointsCreator
  {
  public:
    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoints< Field, dimension > LagrangePoints;
    typedef LagrangePoints LocalCoefficients;

    typedef unsigned int Key;

    template <class Topology>
    struct Setup
    {
      template <int pdim>
      struct InitCodim
      {
        static const unsigned int codim = dimension-pdim;
        static void apply(const unsigned int order,
                          LagrangePoints &points)
        {
          static const unsigned int topologyId = Topology::id;
          typedef GenericGeometry::ReferenceTopologies<dimension> RefElem;
          const unsigned int size = RefElem::get(topologyId).size(codim);
          for (unsigned int i =  0; i < size; ++i)
          {
            const unsigned int subTopologyId = RefElem::get(topologyId).topologyId(codim,i);
            GenericGeometry::IfTopology<Setup<Topology>::template Init,dimension-codim>::apply(subTopologyId,order,i,points);
          }
        }
      };
      template <class SubTopology>
      struct Init
      {
        static const unsigned int codimension = dimension - SubTopology::dimension;
        typedef GenericGeometry::ReferenceMappings< Field, dimension > RefMappings;
        typedef typename RefMappings::Container RefMappingsContainer;
        typedef typename RefMappingsContainer::template Codim< codimension >::Mapping Mapping;
        typedef LabattoInnerPoints<Field,SubTopology> InnerPoints;
        static void apply(const unsigned int order,
                          unsigned int subEntity,
                          LagrangePoints &points)
        {
          unsigned int oldSize = points.size();
          unsigned int size = InnerPoints::size(order);
          points.resize(oldSize+size);
          /*
             std::cout << Topology::name() << " " << SubTopology::name() << " : "
                    << " ( " << codimension <<  " , " << subEntity << " ) "
                    << oldSize << " " << size << std::endl;
           */
          LagrangePoint<Field,dimension> *p = &(points.points_[oldSize]);
          InnerPoints::template setup<dim>(LabattoPoints1D<Field>(order).points_,p);

          const RefMappingsContainer &refMappings = RefMappings::container( Topology::id );
          const Mapping &mapping = refMappings.template mapping< codimension >( subEntity );

          unsigned int nr = 0;
          for ( ; p != &(points.points_[oldSize+size]); ++p )
          {
            /*
               std::cout << "   " << nr << " : "
                      << p->point_ << " -> "
                      << mapping.global( reinterpret_cast< const Dune::FieldVector<Field,dimension-codimension>& >(p->point_) )
                      << std::endl;
             */
            p->point_ = mapping.global( reinterpret_cast< const Dune::FieldVector<Field,dimension-codimension>& >(p->point_) );
            p->localKey_ = LocalKey( subEntity, codimension, nr++ );
          }
        }
      };
    };

    template< class Topology >
    static const LagrangePoints &lagrangePoints ( const Key &order )
    {
      LagrangePoints *lagrangePoints = new LagrangePoints( order, 0 );
      GenericGeometry::ForLoop<Setup<Topology>::template InitCodim,0,dimension>::apply(order,*lagrangePoints);
      return *lagrangePoints;
    }

    template< class T >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      return lagrangePoints< T >( order );
    }

    static void release ( const LagrangePoints &lagrangePoints )
    {
      delete &lagrangePoints;
    }
  };
}
