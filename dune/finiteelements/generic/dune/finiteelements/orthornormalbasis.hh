// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <vector>
namespace Dune
{
  template <class Topology,class F> 
  class OrthonormalBasis 
  {
    typedef OrthonormalBasis<Topology,F> This;
    typedef StandardMonomialBasis<Topology::dimension,F> Basis;

  public:
    typedef typename F Field;

    typedef typename Basis::DomainVector DomainVector;
    typedef typename Basis::RangeVector RangeVector;

    PolynomialBasis () :
    : basis_(), evalBasis_(0)  
    {}

    const int size (unsigned int order) const
    {
      return basis_.size(order);
    }

    void evaluate ( unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basisEval_.resize(size(order));
      basis_.evaluate(order,x,basisEval_);
      mult(basisEval_,values);
    }
  private:
    void mult(const std::vector< RangeVector > &x,
              std::vector< RangeVector > &y) const { 
      // ....
    }
    const Basis basis_;
    mutable std::vector<RangeVector> basisEval_;
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH

