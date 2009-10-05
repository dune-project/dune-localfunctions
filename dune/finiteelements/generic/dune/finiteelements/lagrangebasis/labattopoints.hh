// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

  template <class Field>
  struct LabattoPoints1D
  {
    static const unsigned int size = order-1;
    LabattoPoints1D(unsigned int order) :
    {  
      points_.resize(order-1);
      for (int i=1;i<order;++i)
        points_[i-1] = Field(i)/Field(order);
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
  struct LabattoInnerPoints<Field,Point>
  {
    static const unsigned int dimension = Base::dimension+1;
    template <unsigned int dim>
    static unsigned int size(const unsigned int order)
    {
      return 1;
    }
    template <int dim>
    static void setup(const unsigned int order,
                      const std::vector<double> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      points->point_[0] = 0;
    }
  };
  template <class Field,class Base>
  struct LabattoInnerPoints<Field,Prism<Base> >
  {
    typedef LabattoInnerPoints<Field,Base> LabattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      return LabattoBase::size(order)*(order-1);
    }
    template <unsigned int dim>
    static void setup(const unsigned int order,
                      const std::vector<double> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      assert(dim>=dimension);
      assert(points1D.size()==order-1);
      LabattoBase::setup(order,points1D,points);
      const unsigned int baseSize = LabattoBase::sizei(order);
      for (int q=0;q<order-1;q++)
      {
        for (int i=0;i<baseSize;++i)
        {
          const unsigned int pos = q*baseSize+i;
          for (int d=0;d<dimension-1;++d)
            points[pos].point_[d] = points[i][d];
          points[pos].point_[dimension-1]=points1D[q];
        }
      }
    }
  };
  template <class Field,class Base>
  struct LabattoInnerPoints<Field,Pyramid<Base> >
  {
    typedef LabattoInnerPoints<Field,Base> LabattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      return order*(order-1)/2;
    }
    template <unsigned int dim>
    static void setup(const unsigned int order,
                      const std::vector<double> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      assert(dimension==2);
      assert(points1D.size()==order-1);
      assert(points.size()>size());
      const unsigned int baseSize = LabattoBase::sizei(order);
      unsigned int pos;
      for (int i=2;i<=order;q++)
      {
        for (int j=2;j<=order+1-i;++i,++pos)
        {
          const unsigned int k = order+3-i-j;
          points[pos][0] = 1./3.*(1.+2.*points1D[i]-points1D[j]-points1D[k]);
          points[pos][1] = 1./3.*(1.+2.*points1D[j]-points1D[i]-points1D[k]);
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

    template < class T >
    class SubTopology
    {
      template <int pdim>
      struct Init
      {
        static const unsigned int codim = dimension - pdim;
        static apply(const unsigned int order, std::vector<LagrangePoint> &points) 
        {
          const unsigned int size = GenericGeometry::Size< T, codim>;
          const unsigned int subTopologyId =
                          ReferenceTopologies<dimension>::get(topologyId).topologyId(codim,
        }
      };
    };

    template< class T >
    static const LagrangePoints &lagrangePoints ( const Key &order )
    {
      static const unsigned int topologyId = Topology::id:

      LagrangePoints *lagrangePoints = new LagrangePoints( order, 0 );
      std::vector< LagrangePoint > points = lagrangePoints->points_;
      for (unsigned int codim = 0; codim<=dimension;++codim) {

      }
      GenericGeometry::ForLoop< SubTopology< T >::template Init, 0, dimension >::apply( order, points );
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
