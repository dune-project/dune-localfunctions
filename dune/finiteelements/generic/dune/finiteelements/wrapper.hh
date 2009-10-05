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
struct SimplexWrapper {
  SimplexWrapper(double* storage,int p,
                 int *N)
    : storage_(storage),
      N_(N),
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
  SimplexWrapper(SimplexWrapper& pw,
                 int* N)
    : storage_(pw.storage_),
      N_(N),
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
  ~SimplexWrapper()
  {
    delete [] incr_;
    delete [] row_;
  }
  int p() {
    return p_;
  }
  double* operator()(int r,int c) {
    assert(r<=p_);
    assert(c<=r);
    return row_[r]+incr_[c];
  }
  void fill(double z) {
    for (int r=1; r<=p_; ++r) {
      double *pos0=(*this)(r,0);
      double *pos=pos0;
      for (int c=0; c<r; ++c) {
        for (int i=0; i<N_[c]; ++i,++pos0,++pos) {
          (*pos) = z*(*pos0);
        }
      }
    }
  }
private:
  double* storage_;
  int *N_;
  int *incr_;
  double **row_;  // starting point for Ai0
  int p_;
};
// ********************************************
template <int dim>
struct Simplex {
  static void eval(int p,double *x,double* ret) {
    int N[mat.p()];
    // compute sizes N[k] = |Phi_k|
    SimplexWrapper myMat(ret,p,N);
    Simplex<dim-1>::eval(x,myMat);
    myMat.fill(x[dim]);
  }
  static void eval(double *x,SimplexWrapper& mat) {
    int N[mat.p()];
    // compute sizes N[k] = |Phi_k|
    SimplexWrapper myMat(mat,N);
    Simplex<dim-1>::eval(x,myMat);
    myMat.fill(x[dim]);
  }
};
