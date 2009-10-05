// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <cassert>
/**
 * Use a given storage to define structure of the form
 *   A00
 *   A10 A11
 *   A20 A21 A22
 *   ...
 *   Ap0 Ap1 ... App
 * Here A(r,c) is double[Nc] where N0,..,Np are given in constructor.
 * Assumption is that (Ar0,...,Arr) are consequtive in memory but
 * there might be a gap between Arr und A(r+1)0.
 * Methods on Wrapper:
 * double* operator()(int p,int i) -> Api
 *
 **/
struct PyramidWrapper {
  PyramidWrapper(double* storage,
                 int* N,int p)
    : storage_(storage),
      row_(new double*[p+1]),
      incr_(new int[p+1]), // Arc = row[p]+incr[c]
      p_(p)
  {
    incr_[0]=0;
    row_[0]=storage_;
    for (int r=0; r<p_; ++r) {
      incr_[r+1]=incr_[r]+N[r];
      row_[r+1]=row_[r]+incr_[r+1];
    }
  }
  PyramidWrapper(PyramidWrapper& pw,
                 int* N)
    : storage_(pw.storage_),
      row_(new double*[pw.p_+1]),      // are Arr from pw
      incr_(new int[pw.p_+1]),     // increments are new
      p_(pw.p_)
  {
    incr_[0]=0;
    row_[0]=pw(0,0);
    for (int r=0; r<p_; ++r) {
      incr_[r+1]=incr_[r]+N[r];
      row_[r+1]=pw(r+1,r+1);
    }
  }
  ~PyramidWrapper()
  {
    delete [] incr_;
    delete [] row_;
  }
  double* operator()(int r,int c) {
    assert(r<=p);
    assert(c<=r);
    return row_[r]+incr_[c];
  }
private:
  double* storage_;
  int *incr_;
  double **row_;  // starting point for Ai0
  int p_;
};
