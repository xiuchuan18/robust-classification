// tpsc.cpp  -- Rcpp core for the S-RoLLR classifier
#include <Rcpp.h>
using namespace Rcpp;

static const double LOG_1e100 = -230.2585092994046;   // log(1e-100)
static const double LOG_2PI   = 1.8378770664093453;   // log(2*pi)

// Precomputed, delta-dependent pieces of the standard Student-t log-density,
// so the (expensive) lgamma calls are evaluated once per feature, not per point.
struct TLogConst {
  bool   normal;   // delta >= 1000 then standard normal base
  double cst;      // lgamma((d+1)/2) - lgamma(d/2) - 0.5*log(d*pi)
  double a;        // (delta+1)/2
  double delta;
  TLogConst(double d) : delta(d) {
    normal = (d >= 1000.0);
    if (normal) { cst = -0.5 * LOG_2PI; a = 0.0; }
    else {
      cst = R::lgammafn((d + 1.0) / 2.0) - R::lgammafn(d / 2.0)
          - 0.5 * std::log(d * M_PI);
      a   = (d + 1.0) / 2.0;
    }
  }
  inline double logpdf(double z) const {
    if (normal) return cst - 0.5 * z * z;
    return cst - a * std::log1p(z * z / delta);
  }
};

// Precomputed pieces of a single TPSC log-density (fixed theta,w,sigma,delta).
struct TpscConst {
  bool   ok;
  double theta, s1, s2, l2w, l2mw, ls1, ls2;
  TLogConst tc;
  TpscConst(double theta_, double w, double sigma, double delta)
      : tc(delta) {
    ok = (w > 0.0 && w < 1.0 && sigma > 0.0 && delta > 0.0);
    theta = theta_;
    if (ok) {
      s1  = sigma * std::sqrt(w / (1.0 - w));
      s2  = sigma * std::sqrt((1.0 - w) / w);
      l2w = std::log(2.0 * w);      ls1 = std::log(s1);
      l2mw= std::log(2.0 * (1.0-w)); ls2 = std::log(s2);
    }
  }
  inline double logpdf(double x) const {
    if (!ok) return LOG_1e100;
    double lf = (x < theta) ? (l2w  - ls1 + tc.logpdf((x - theta) / s1))
                            : (l2mw - ls2 + tc.logpdf((x - theta) / s2));
    return (lf < LOG_1e100) ? LOG_1e100 : lf;
  }
};

// Self-contained scalar version (for varying params), kept for d_tpsc_cpp.
inline double tpsc_logpdf1(double x, double theta, double w,
                           double sigma, double delta) {
  return TpscConst(theta, w, sigma, delta).logpdf(x);
}

// [[Rcpp::export]]
NumericVector d_tpsc_cpp(NumericVector x, double w, double theta,
                         double sigma, double delta) {
  R_xlen_t n = x.size();
  NumericVector out(n);
  if (w <= 0.0 || w >= 1.0 || sigma <= 0.0 || delta <= 0.0) return out; // zeros
  for (R_xlen_t i = 0; i < n; ++i)
    out[i] = std::exp(tpsc_logpdf1(x[i], theta, w, sigma, delta));
  return out;
}

// Negative log-likelihood, parameter order (theta, w, sigma, delta) to match optim().
// [[Rcpp::export]]
double tpsc_nll_cpp(NumericVector params, NumericVector data) {
  if (params.size() != 4) return R_PosInf;
  double theta = params[0], w = params[1], sigma = params[2], delta = params[3];
  if (!R_finite(theta) || !R_finite(w) || !R_finite(sigma) || !R_finite(delta))
    return R_PosInf;
  if (w <= 1e-6 || w >= 1.0 - 1e-6 || sigma <= 1e-8 || delta <= 0.1 || delta >= 1000.0)
    return R_PosInf;
  TpscConst tp(theta, w, sigma, delta);   // lgamma etc. computed once
  double s = 0.0;
  R_xlen_t n = data.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = data[i];
    if (!R_finite(xi)) continue;
    s += tp.logpdf(xi);
  }
  if (!R_finite(s)) return R_PosInf;
  return -s;
}

// T_ij = log f_{1j}(X_ij) - log f_{0j}(X_ij).
// P0, P1 are p*4 matrices of per-feature params (theta, w, sigma, delta)
// for class 0 and class 1 respectively. Vectorized in C++.
// [[Rcpp::export]]
NumericMatrix tpsc_transform_cpp(NumericMatrix X, NumericMatrix P0, NumericMatrix P1) {
  R_xlen_t n = X.nrow(), p = X.ncol();
  NumericMatrix Tr(n, p);
  for (R_xlen_t j = 0; j < p; ++j) {
    TpscConst c0(P0(j,0), P0(j,1), P0(j,2), P0(j,3));   // constants once per feature
    TpscConst c1(P1(j,0), P1(j,1), P1(j,2), P1(j,3));
    for (R_xlen_t i = 0; i < n; ++i) {
      double x = X(i,j);
      Tr(i,j) = c1.logpdf(x) - c0.logpdf(x);
    }
  }
  return Tr;
}

// Total log-likelihood of each row of X under one class's per-feature params
// [[Rcpp::export]]
NumericVector tpsc_rowloglik_cpp(NumericMatrix X, NumericMatrix P) {
  R_xlen_t n = X.nrow(), p = X.ncol();
  NumericVector out(n);
  std::vector<TpscConst> cc;               // constants once per feature
  cc.reserve(p);
  for (R_xlen_t j = 0; j < p; ++j)
    cc.emplace_back(P(j,0), P(j,1), P(j,2), P(j,3));
  for (R_xlen_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (R_xlen_t j = 0; j < p; ++j) s += cc[j].logpdf(X(i,j));
    out[i] = s;
  }
  return out;
}
