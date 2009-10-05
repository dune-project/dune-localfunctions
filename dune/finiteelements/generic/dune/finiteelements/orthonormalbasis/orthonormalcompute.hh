// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <cassert>
#include <alglib/qr.h>
#include <alglib/sevd.h>
#include <dune/fem/space/lagrangespace/genericlagrangepoints.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
template <int d>
struct MultiIndex { // FieldVector<int,d> {
  int absval;
  std::vector<int> dat;
  MultiIndex(int v=0) : dat(d) {
    set(v);
  }
  MultiIndex& operator=(int v) {
    set(v);
    return *this;
  }
  MultiIndex(const MultiIndex<d>& other) :
    absval(other.absval),
    dat(other.dat) {
    assert(dat.size()==other.dat.size());
  }
  MultiIndex& operator=(const MultiIndex& other) {
    dat=other.dat;
    absval = other.absval;
    return *this;
  }
  void set(int v) {
    for (int i=0; i<d; i++) {
      set(i,v);
    }
    absval = v*d;
  }
  void set(int i,int v) {
    absval += v - dat[i];
    dat[i] = v;
  }
  const int& operator[](int i) const {
    return dat[i];
  }
  bool incr(int ord) {
    for (int i=0; i<d; i++)
      if (dat[i]<ord && absval<ord) {
        dat[i]++;
        absval++;
        return true;
      } else {
        absval -= dat[i];
        dat[i]=0;
      }
    return false;
  }
  int abs() {
    int ret = dat[0];
    for (int i=1; i<d; i++)
      ret += dat[i];
    assert(absval==ret);
    return absval;
  }
  int abs(int n) {
    int ret = dat[0];
    for (int i=1; i<n; i++)
      ret += dat[i];
    return ret;
  }
  inline MultiIndex<d>& operator+=(const MultiIndex<d>& m) {
    for (int i=0; i<d; i++)
      dat[i] += m.dat[i];
    absval += m.absval;
    return *this;
  }
};
template<int d>
inline bool operator<(const MultiIndex<d>& x, const MultiIndex<d>& y) {
  return (x.dat < y.dat);
}
template<int n>
std::ostream& operator<< (std::ostream& s, const MultiIndex<n>& v)
{
  for (int i=0; i<n; i++)
    s << ((i>0) ? " " : "") << v[i];
  return s;
}
template <int d,int ord>
struct Monoms {
  typedef enum {lexi=1,absvalue=2} order_t;
  typedef Dune::GenericLagrangePoint<
      typename Dune::GeometryWrapper<Dune::GeometryType :: simplex,d>::GenericGeometryType,ord>
  PointType;
  enum {size=PointType::numLagrangePoints};
  static void multiindex(int idx,MultiIndex<d>& coord,
                         order_t order=absvalue) {
    if (order==absvalue)
      idx = size-1-idx;
    Dune::FieldVector<double,d> dcoord;
    PointType points(idx);
    points.local(dcoord);
    dcoord *= double(ord);
    if (order==lexi) {
      for (int i=0; i<d; i++)
        coord.set(i,int(dcoord[d-1-i]));
    }
    else {
      double bary = dcoord[d-1];
      for (int i=d-2; i>=0; --i) {
        bary += dcoord[i];
        coord.set(i+1,int(dcoord[i]));
      }
      coord.set(0,ord-int(bary));
    }
  }
  Monoms() {
    for (int i=0; i<size; i++)
      values_[i] = -1.;
  }
  double values_[size];
};

typedef amp::ampf<Precision> scalar_t;
typedef ap::template_1d_array< scalar_t > vec_t;
typedef ap::template_2d_array< scalar_t > mat_t;

scalar_t factorial(int start,int end) {
  scalar_t ret(1);
  for (int j=start; j<=end; j++)
    ret*=j;
  return ret;
}

template <int d,int dim>
struct IntegralHelper {
  static void compute(MultiIndex<dim>& geo,
                      MultiIndex<dim>& alpha,scalar_t& p,scalar_t& q) {
    IntegralHelper<d-1,dim>::compute(geo,alpha,p,q);
    if (geo[d-1] == 1) {
      p *= factorial(1,alpha[d-1]);
      q *= factorial(d+alpha.abs(d-1),d+alpha.abs(d-1)+alpha[d-1]);
    } else if (geo[d-1] == 2) {
      p *= 1.;
      q *= alpha[d-1]+1.;
    }
    else abort();
  }
};
template <int dim>
struct IntegralHelper<1,dim> {
  static void compute(MultiIndex<dim>& geo,
                      MultiIndex<dim>& alpha,scalar_t& p,scalar_t& q) {
    if (geo[0] == 1) {
      p = 1.;
      q = alpha[0]+1;
    } else if (geo[0] == 2) {
      p = 1.;
      q = alpha[0]+1;
    }
    else abort();
  }
};
template <int dim,int ord>
struct Integral {
  static void compute(MultiIndex<dim>& geo,
                      int i,int j,scalar_t& p,scalar_t& q) {
    MultiIndex<dim> coordi,coordj;
    Monoms<dim,ord>::multiindex(i,coordi);
    Monoms<dim,ord>::multiindex(j,coordj);
    coordi += coordj;
    IntegralHelper<dim,dim>::compute(geo,coordi,p,q);
  }
};
/**************************************************/
template <int dim,int ord>
struct CalcCoeffsBase {
  enum {N = Monoms<dim,ord>::size};
  vec_t d;
  mat_t S;
  mat_t& res;
  CalcCoeffsBase(mat_t& pres) : res(pres) {
    res.setbounds(1,N,1,N);
    S.setbounds(1,N,1,N);
    d.setbounds(1,N);
  }
  virtual ~CalcCoeffsBase() {}
  void test() {
    // Nun schauen wir noch, ob wirklich C^T S C = E gilt
    amp::ampf<Precision> prod;
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        prod = 0;
        for (int k=1; k<=N; k++) {
          for (int l=1; l<=N; l++) {
            prod += res(l,i)*S(l,k)*res(k,j);
          }
        }
        assert((i==j) ? 1 : // fabs(prod.toDouble()-1)<1e-10 :
               fabs(prod.toDouble())<1e-10);
        // std::cout << "(" << i << "," << j << ")" << prod.toDouble() << std::endl;
      }
    }
  }
  virtual void compute() = 0;
  void compute(MultiIndex<dim>& geo) {
    // Aufstellen der Matrix fuer die Bilinearform xSy
    // S_ij = int_A x^(i+j)
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        amp::ampf<Precision> p,q;
        Integral<dim,ord>::compute(geo,i-1,j-1,p,q);
        S(i,j) = p;
        S(i,j) /= q;
      }
    }
    compute();
    test();
    return;
  }
};
template <int dim,int ord>
struct CalcCoeffsGM : public CalcCoeffsBase<dim,ord> {
  enum {N = Monoms<dim,ord>::size};
  CalcCoeffsGM(mat_t& pres) : CalcCoeffsBase<dim,ord>(pres) {}
  virtual ~CalcCoeffsGM() {}
  void sprod(int col1,int col2,
             scalar_t& ret) {
    ret = 0;
    for (int k=1; k<=col1; k++) {
      for (int l=1; l<=col2; l++) {
        ret += this->res(l,col2)*this->S(l,k)*this->res(k,col1);
      }
    }
  }
  void compute() {
    scalar_t s;
    for (int i=1; i<=N; ++i)
      for (int j=1; j<=N; ++j)
        this->res(i,j) = (i==j) ? 1 : 0;
    sprod(1,1,s);
    s = 1./amp::sqrt<Precision>(s);
    vmul(this->res.getcolumn(1,1,1),s);
    for (int i=1; i<N; i++) {
      for (int k=1; k<=i; ++k) {
        sprod(i+1,k,s);
        vsub(this->res.getcolumn(i+1,1,i+1),
             this->res.getcolumn(k,  1,i+1),s);
      }
      sprod(i+1,i+1,s);
      s = 1./amp::sqrt<Precision>(s);
      vmul(this->res.getcolumn(i+1,1,i+1),s);
    }
  }
};
template <int dim,int ord>
struct CalcCoeffsQR : public CalcCoeffsBase<dim,ord> {
  enum {N = Monoms<dim,ord>::size};
  mat_t a,q,r,z,Ssqrt,Ssqrtinv,tmp,tmpinv;
  CalcCoeffsQR(mat_t& pres) : CalcCoeffsBase<dim,ord>(pres) {
    q.setbounds(1,N,1,N);
    r.setbounds(1,N,1,N);
    z.setbounds(1,N,1,N);
    Ssqrt.setbounds(1,N,1,N);
    tmp.setbounds(1,N,1,N);
    Ssqrtinv.setbounds(1,N,1,N);
    tmpinv.setbounds(1,N,1,N);
  }
  virtual ~CalcCoeffsQR() {}
  void compute() {

    /*
       for (int i=1;i<=N;++i) {
       for (int j=1;j<=N;j++)
        std::cout << S(i,j).toDouble() << "\t\t" << std::flush;
       std::cout << std::endl;
       }
     */

    // Berechne die Eigenwerte d_1,..,d_N und ein ONS von
    // Rechtseigenvektoren
    sevd::symmetricevd(this->S,N,1,true,this->d,z);
    // Berechne SQRT(S) und SQRT(S)^(-1)
    // Dazu: S=R^T D R und S^(-1) = R^T D^(-1) R
    for (int i=1; i<=N; ++i) {
      this->d(i) = amp::sqrt<Precision>(this->d(i));
    }
    for (int i=1; i<=N; ++i)
      for (int j=1; j<=N; j++) {
        Ssqrt(i,j) = 0;
        Ssqrtinv(i,j) = 0;
        tmp(i,j) = 0;
        tmpinv(i,j) = 0;
      }
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        tmp(i,j) = this->d(i)*z(j,i);
        tmpinv(i,j) = z(j,i)/this->d(i);
      }
    }
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        for (int k=1; k<=N; k++) {
          Ssqrt(i,j) = Ssqrt(i,j) + z(i,k)*tmp(k,j);
          Ssqrtinv(i,j) = Ssqrtinv(i,j) + z(i,k)*tmpinv(k,j);
        }
      }
    }
    // Testen das alles so weit stimmt....
    for (int i=1; i<=N; ++i)
      for (int j=1; j<=N; j++) {
        tmp(i,j) = 0;
        tmpinv(i,j) = 0;
      }
    // SQRT(S), SQRT(S)^(-1) sind symmetrisch ...
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        assert(fabs(Ssqrt(i,j).toDouble() - Ssqrt(j,i).toDouble())<1e-10);
        assert(fabs(Ssqrtinv(i,j).toDouble() - Ssqrtinv(j,i).toDouble())<1e-10);
        for (int k=1; k<=N; k++) {
          tmp(i,j) = tmp(i,j) + Ssqrt(k,i)*Ssqrt(k,j);
          tmpinv(i,j) = tmpinv(i,j) + Ssqrt(k,i)*Ssqrtinv(k,j);
        }
      }
    }
    // ... und SQRT(S)^2 = S
    // ... und SQRT(S)*SQRT(S)^(-1) = E
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        assert(fabs(tmp(i,j).toDouble() - this->S(i,j).toDouble())<1e-10);
        assert((i==j) ? (tmpinv(i,j).toDouble()-1)<1e-10 : (tmpinv(i,j).toDouble())<1e-10);
      }
    }
    // Nun berechne die QR Zerlegung von SQRT(S)
    qr::qrdecompositionunpacked(Ssqrt,N,N,q,r);
    // Nun die gesuchten Koeffizienten:
    // C = SQRT(S)^(-1)q
    // Es gilt jetzt:
    // SQRT(S) = qr
    // ->   E = SQRT(S)^(-1)qr = Cr
    // Da: q^T = q^(-1) bzw. q^T SQRT(S) = r
    // ->   E = Cr = Cq^T SQRT(S)
    // C^T S C = q^T SQRT(S)^(-T)S SQRT(S)^(-1) q = q^T S^(-1)S q = q^Tq = E
    // d.h. die Spalten sind orthogonal bzw. der Bilinearform S!
    for (int i=1; i<=N; ++i) {
      for (int j=1; j<=N; j++) {
        this->res(i,j) = 0;
        for (int k=1; k<=N; k++) {
          this->res(i,j) = this->res(i,j) + Ssqrtinv(k,i)*q(k,j);
        }
      }
    }
    // ... wir wollen noch einen positiven f\"uhrenden Koeffizienten ...
    for (int i=1; i<=N; ++i) {
      if (this->res(i,i).isNegativeNumber())
        for (int j=1; j<=N; j++) {
          this->res(j,i) *= -1.0;
        }
      // ... event. Normieren ...
      /*
         for (int j=1;j<=N;j++) {
         res(j,i) = res(j,i) / res(i,i);
         }
       */
    }
  }
};
template <int dim,int ord>
struct CalcCoeffs {
  enum {N = Monoms<dim,ord>::size};
  mat_t res;
  CalcCoeffsBase<dim,ord> *method;
  CalcCoeffs(int use_method) {
    if (use_method==1)
      method = new CalcCoeffsGM<dim,ord>(res);
    else
      method = new CalcCoeffsQR<dim,ord>(res);
  }
  ~CalcCoeffs() {
    delete method;
  }
  void compute(MultiIndex<dim>& geo) {
    method->compute(geo);
  }
};
