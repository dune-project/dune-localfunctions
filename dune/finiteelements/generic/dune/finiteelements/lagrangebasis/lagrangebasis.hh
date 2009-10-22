// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/finiteelements/common/matrix.hh>
#include <dune/finiteelements/generic/basismatrix.hh>

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

#if 0
  template< template <class,unsigned int> class LP,
      class Topology, class scalar_t >
  struct LagrangeMatrix
  {
    static const unsigned int dimension = Topology::dimension;

    typedef LFEMatrix< scalar_t > mat_t;
    typedef LagrangeInterpolationFactory< LP,dimension,scalar_t > LocalInterpolationFactory;
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

    unsigned int cols () const
    {
      return matrix_.cols();
    }

    unsigned int rows () const
    {
      return matrix_.rows();
    }

    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      for (unsigned int i=0; i<cols(); ++i)
        field_cast(matrix_(i,row),vec[i]);
    }

    /*
       scalar_t& operator() ( int r, int c ) const
       {
       return Dune::FieldMatrix< scalar_t, 1, 1 >( matrix_( c, r ) );
       }
     */

    void print ( std::ostream &out, const unsigned int N = rows() ) const
    {
      for( unsigned int i = 0; i < N; ++i )
      {
        out << "Polynomial : " << i << std::endl;
        for( unsigned int j = 0; j <cols(); ++j )
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
#endif


  // LagrangeBasisFactory
  // --------------------

  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisFactory;
  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisTraits
  {
    typedef Dune::MonomialBasisProvider< dim, SF > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Object MonomialBasis;
    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< SF, 1 > > Basis;

    typedef Dune::MonomialBasisProvider< dim, CF > PreBasisFactory;
    typedef typename PreBasisFactory::Object PreBasis;
    typedef LagrangeInterpolationFactory< LP, dim, CF > InterpolationFactory;
    typedef typename InterpolationFactory::Object Interpolation;

    static const unsigned int dimension = dim;
    typedef const Basis Object;
    typedef unsigned int Key;
    typedef LagrangeBasisFactory<LP,dim,SF,CF> Factory;
  };

  template< template <class,unsigned int> class LP,
      unsigned int dim, class SF, class CF >
  struct LagrangeBasisFactory
    : public TopologyFactory< LagrangeBasisTraits< LP,dim,SF,CF > >
  {
    typedef LagrangeBasisTraits<LP,dim,SF,CF> Traits;
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

      const typename Traits::PreBasis *preBasis = Traits::PreBasisFactory::template create<Topology>( order );
      const typename Traits::Interpolation *interpol = Traits::InterpolationFactory::template create<Topology>( order );
      BasisMatrix< typename Traits::PreBasis,
          typename Traits::Interpolation,
          ComputeField > matrix( *preBasis, *interpol );
      basis->fill( matrix );
      Traits::PreBasisFactory::release(preBasis);
      Traits::InterpolationFactory::release(interpol);

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
