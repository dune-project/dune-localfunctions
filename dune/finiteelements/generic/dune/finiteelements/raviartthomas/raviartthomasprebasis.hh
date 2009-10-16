// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASPREBASIS_HH
#define DUNE_RAVIARTTHOMASPREBASIS_HH
#include <fstream>
#include <utility>

#include <dune/common/forloop.hh>

#include <dune/finiteelements/common/matrix.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>

#include <dune/finiteelements/common/localcoefficients.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#if HAVE_ALGLIB
#include <dune/finiteelements/lagrangebasis/lobattopoints.hh>
#endif
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/quadrature/subquadrature.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>
namespace Dune
{
  template <class Topology, class Field>
  struct RTPreBasisCreator
  {
    enum {dim = Topology::dimension};
    typedef MonomialBasis<Topology,Field> MBasis;
    typedef StandardEvaluator<MBasis> EvalMBasis;
    typedef PolynomialBasisWithMatrix<EvalMBasis,SparseCoeffMatrix<Field,dim> > Basis;
    struct RTVecMatrix
    {
      typedef MultiIndex<dim> MI;
      typedef MonomialBasis<Topology,MI> MIBasis;
      RTVecMatrix(unsigned int order)
      {
        MIBasis basis(order+1);
        FieldVector< MI, dim > x;
        for( int i = 0; i < dim; ++i )
          x[ i ].set( i, 1 );
        std::vector< MI > val( basis.size() );
        basis.evaluate( x, val );

        col_ = basis.size();
        unsigned int notHomogen = 0;
        if (order>0)
          notHomogen = basis.sizes()[order-1];
        unsigned int homogen = basis.sizes()[order]-notHomogen;
        row_ = (notHomogen*dim+homogen*(dim+1))*dim;
        row1_ = basis.sizes()[order]*dim*dim;
        mat_ = new Field*[row_];
        int row = 0;
        for (unsigned int i=0; i<notHomogen+homogen; ++i)
        {
          for (unsigned int r=0; r<dim; ++r)
          {
            for (unsigned int rr=0; rr<dim; ++rr)
            {
              mat_[row] = new Field[col_];
              for (unsigned int j=0; j<col_; ++j)
              {
                mat_[row][j] = 0.;
              }
              if (r==rr)
                mat_[row][i] = 1.;
              ++row;
            }
          }
        }
        for (unsigned int i=0; i<homogen; ++i)
        {
          for (unsigned int r=0; r<dim; ++r)
          {
            mat_[row] = new Field[col_];
            for (unsigned int j=0; j<col_; ++j)
            {
              mat_[row][j] = 0.;
            }
            unsigned int w;
            MI xval = val[notHomogen+i];
            xval *= x[r];
            for (w=homogen+notHomogen; w<val.size(); ++w)
            {
              if (val[w] == xval)
              {
                mat_[row][w] = 1.;
                break;
              }
            }
            assert(w<val.size());
            ++row;
          }
        }
      }
      ~RTVecMatrix()
      {
        for (unsigned int i=0; i<rowSize(); ++i) {
          delete [] mat_[i];
        }
        delete [] mat_;
      }
      unsigned int colSize(int row) const {
        return col_;
      }
      unsigned int rowSize() const {
        return row_;
      }
      const Field operator() ( int r, int c ) const
      {
        assert(r<(int)row_ && c<(int)col_);
        return mat_[r][c];
      }
      unsigned int row_,col_,row1_;
      Field **mat_;
    };
    static Basis &basis(unsigned int order)
    {
      RTVecMatrix vecMatrix(order);
      MBasis basis(order+1);
      Basis *tmBasis = new Basis(basis);
      tmBasis.fill(vecMatrix);
      return *tmBasis;
    }
    static void release ( const Basis &basis )
    {
      delete &basis;
    }
  };
  // -----------------------------------------
  // Basis
  // -----------------------------------------
  template <class Field,unsigned int dim>
  struct RaviartThomasInitialBasis
  {
    static const unsigned int dimension = dim;
    // typedef MonomialBasisProvider<dimension,Field> TestBasisProvider;
    // typedef MonomialBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    typedef OrthonormalBasisProvider<dimension,Field> TestBasisProvider;
    typedef OrthonormalBasisProvider<dimension-1,Field> TestFaceBasisProvider;
    // typedef LagrangeBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LagrangeBasisProvider<dimension-1,Field> TestFaceBasisProvider;
#if HAVE_ALGLIB
    // typedef LobattoBasisProvider<dimension,Field> TestBasisProvider;
    // typedef LobattoBasisProvider<dimension-1,Field> TestFaceBasisProvider;
#endif

    typedef typename TestBasisProvider::Basis TestBasis;
    typedef typename TestFaceBasisProvider::Basis TestFaceBasis;

    typedef MonomialBasisProvider<dimension,Field> MBasisProvider;
    typedef typename MBasisProvider::Basis MBasis;

    struct FaceStructure
    {
      FaceStructure( const TestFaceBasis &tfb, const Dune::FieldVector<Field,dimension> &n )
        : basis(tfb), normal(n)
      {}
      const TestFaceBasis &basis;
      const Dune::FieldVector<Field,dimension> &normal;
      template < class Topology >
      struct Creator
      {
        template < int face >
        struct GetCodim
        {
          typedef typename GenericGeometry::SubTopology<Topology,1,face>::type FaceTopology;
          static void apply( const unsigned int order,
                             std::vector<FaceStructure*> &faceStructure )
          {
            faceStructure.push_back( new FaceStructure (
                                       TestFaceBasisProvider::template basis<FaceTopology>(order),
                                       GenericGeometry::ReferenceElement<Topology,Field>::integrationOuterNormal(face) ) );
          }
        };
      };
    };
    template < class Topology >
    struct Creator
    {
      typedef RTPreBasisCreator<Topology,Field> PreBasisCreator;
      typedef typename PreBasisCreator::Basis PreBasis;

      static const TestBasis &testBasis(unsigned int order)
      {
        return TestBasisProvider::template basis<Topology>(order-1);
      }
      static const MBasis &mBasis(unsigned int order)
      {
        return MBasisProvider::template basis<Topology>(order+1);
      }
      static const PreBasis &preBasis(unsigned int order)
      {
        return PreBasisCreator::basis(order);
      }
    };
  };
}
#endif // DUNE_RAVIARTTHOMASPREBASIS_HH
