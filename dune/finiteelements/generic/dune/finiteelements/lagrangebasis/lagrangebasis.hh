// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/finiteelements/common/matrix.hh>

#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/generic/multiindex.hh>
#include <dune/finiteelements/generic/monomialbasis.hh>
#include <dune/finiteelements/generic/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>

namespace Dune
{

  // LagrangeMatrix
  // --------------

  template< class Topology, class scalar_t,
      class LPFactory >
  struct LagrangeMatrix
  {
    static const unsigned int dimension = Topology::dimension;

    typedef LPFactory LagrangeCoefficientsFactory;
    typedef LFEMatrix< scalar_t > mat_t;
    typedef LagrangeInterpolationFactory< LagrangeCoefficientsFactory > LocalInterpolationFactory;
    typedef typename LocalInterpolationFactory::Object LocalInterpolation;

    explicit LagrangeMatrix( const unsigned int order )
    {
      Dune::MonomialBasis< Topology, scalar_t > basis( order );
      const LocalInterpolation *localInterpolation
        = LocalInterpolationFactory::template create< Topology >( order );
      localInterpolation->interpolate( basis, matrix_ );
      LocalInterpolationFactory::release( localInterpolation );
      if ( !matrix_.invert() )
      {
        DUNE_THROW(MathError, "While computing LagrangeBasis a singular matrix was constructed!");
      }
    }

    unsigned int colSize( const unsigned int row ) const
    {
      return matrix_.cols();
    }

    unsigned int rowSize () const
    {
      return matrix_.rows();
    }

    Dune::FieldMatrix< scalar_t, 1, 1 > operator() ( int r, int c ) const
    {
      return Dune::FieldMatrix< scalar_t, 1, 1 >( matrix_( c, r ) );
    }

    void print ( std::ostream &out, const unsigned int N = rowSize() ) const
    {
      for( unsigned int i = 0; i < N; ++i )
      {
        out << "Polynomial : " << i << std::endl;
        for( unsigned int j = 0; j <colSize( i ); ++j )
        {
          double v = matrix_( j, i ).toDouble();
          if( fabs( v ) < 1e-20 )
            out << 0 << "\t\t" << std::flush;
          else
          {
            scalar_t v = matrix_( j, i );
            out << v << "\t\t" << std::flush;
          }
        }
        out << std::endl;
      }
    }

  private:
    mat_t matrix_;
  };



  // LagrangeBasisFactory
  // --------------------


  template< unsigned int dim,
      template <class Field,unsigned int> class LPFactory,
      class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LagrangeBasisFactory;
  template< unsigned int dim,
      template <class Field,unsigned int> class LPFactory,
      class SF, class CF >
  struct LagrangeBasisTraits
  {
    typedef Dune::MonomialBasisProvider< dim, SF > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Object MonomialBasis;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< SF, 1 > > Basis;

    static const unsigned int dimension = dim;
    typedef const Basis Object;
    typedef unsigned int Key;
    typedef LagrangeBasisFactory<dim,LPFactory,SF,CF> Factory;
  };

  template< unsigned int dim,
      template <class Field,unsigned int> class LPFactory,
      class SF, class CF >
  struct LagrangeBasisFactory
    : public TopologyFactory< LagrangeBasisTraits< dim,LPFactory,SF,CF > >
  {
    typedef LagrangeBasisTraits<dim,LPFactory,SF,CF> Traits;
    static const unsigned int dimension = dim;
    typedef SF StorageField;
    typedef CF ComputeField;
    typedef typename Traits::Basis Basis;

    typedef typename Traits::Object Object;
    typedef typename Traits::Key Key;

    template< class Topology >
    static Object *createObject ( const Key &order )
    {
      const typename Traits::MonomialBasis *monomialBasis = Traits::MonomialBasisProvider::template create< Topology >( order );
      Basis *basis = new Basis( *monomialBasis );
      LagrangeMatrix< Topology, ComputeField, LPFactory<ComputeField,dim> > matrix( order );
      basis->fill( matrix );
#if GLFEM_BASIS_PRINT
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,1> > BasisMI;
        const MBasisMI &_mBasisMI = Dune::MonomialBasisProvider<dimension,MIField>::template create<Topology>(order);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order << ".basis";
        std::ofstream out(name.str().c_str());
        basisPrint<0>(out,basisMI);
      }
#endif
      return basis;
    }
  };
}

#endif // #ifndef DUNE_LAGRANGEBASIS_HH
