// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_COEFFMATRIX_HH
#define DUNE_COEFFMATRIX_HH
#include <cassert>
#include <dune/alglib/libs/ap.h>
#include <dune/alglib/libs/amp.h>

namespace Dune
{
  template <class F>
  class CoeffMatrix {
    typedef F Field;
  public:
    CoeffMatrix()
      : coeff_(0),
        rows_(0),
        numRows_(0)
    {}
    ~CoeffMatrix() {
      delete [] coeff_;
      delete [] rows_;
    }
    const int size () const {
      return numRows_;
    }
    template <class RangeVector>
    void mult(const std::vector< RangeVector > &x,
              std::vector< RangeVector > &y) const {
      size_t numLsg = y.size();
      assert(numLsg<=numRows_);
      Field *row = rows_[0];
      for (int r=0; r<numLsg; ++r) {
        y[r] = Field(0);
        RangeVector *const itx = (&x[0]);
        for (; row != rows_[r+1]; ++row, ++itx) {
          y[r] += (*row)*(*itx);
        }
      }
    }
    template <class FullMatrix>
    void fill(const FullMatrix& mat) {
      numRows_ = mat.rowSize();
      int size = 0;
      for (int r=0; r<numRows_; ++r)
        size += mat.colSize(r);
      coeff_ = new Field[size];
      rows_ = new Field*[numRows_+1];
      rows_[0] = coeff_;
      for (int r=0; r<numRows_; ++r) {
        rows_[r+1] = rows_[r] + mat.colSize(r);
        int c = 0;
        for (Field *it = rows_[r]; it != rows_[r+1]; ++it,++c) {
          mat.set(r,c,*it);
        }
      }
    }
  private:
    Field *coeff_;
    Field **rows_;
    unsigned int numRows_;
  };
}

#endif // DUNE_COEFFMATRIX_HH
