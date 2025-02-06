#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

vec seq_C(const double& a, const double& b, const double& by) {
  int n = round((b - a) / by) + 1;
  vec out = linspace(a, b, n);
  return out;
}

vec PSurfaceH0(const rowvec& point,
               const vec& y_lower,
               const vec& y_upper,
               const double& sigma,
               const double& maxP) {
  double y = point[1];
  vec out2 = normcdf(y_upper, y, sigma) - normcdf(y_lower, y, sigma);
  double wgt = y / maxP * (y < maxP) + (y >= maxP);
  return wgt * out2;
}

vec PSurfaceHk(const rowvec& point,
                     const vec& y_lower,
                     const vec& y_upper,
                     const vec& x_lower,
                     const vec& x_upper,
                     const double& sigma,
                     const double& maxP) {
  double x = point[0];
  double y = point[1];
  vec out1 = normcdf(x_upper, x, sigma) - normcdf(x_lower, x, sigma);
  vec out2 = normcdf(y_upper, y, sigma) - normcdf(y_lower, y, sigma);
  double wgt = y / maxP * (y < maxP) + (y >= maxP);
  return wgt * vectorise(out1 * out2.t()); // compute outer product and flatten
}

// [[Rcpp::export]]
NumericVector computePersistenceImage(const mat& D,
                    const int& homDim,
                    const vec& xSeq,
                    const vec& ySeq,
                    const double& sigma) {

  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  // If there are no matching rows, return a zero vector
  if (indices.n_elem == 0) return NumericVector((xSeq.n_elem - 1)*(ySeq.n_elem - 1));

  int resP = ySeq.n_elem - 1;
  int resB = xSeq.n_elem - 1;
  int n_rows = indices.n_elem;

  if (n_rows == 0) return (resB > 0) ? NumericVector(resB * resP) : NumericVector(resP);

  mat D_=D.submat(indices,uvec{1,2});

  double minP = ySeq[0];
  double maxP = ySeq[ySeq.n_elem - 1];
  double dy = (maxP - minP) / resP;
  vec y_lower = seq_C(minP, maxP - dy, dy);
  vec y_upper = y_lower + dy;

  int Nsize;
  double sumB = sum(abs(diff(D_.col(0))));
  if ((homDim == 0) && (sumB == 0)) {
    Nsize = resP;
  } else {
    Nsize = resB * resP;
  }

  mat Psurf_mat(Nsize, n_rows);

  if (Nsize == resP) {
    for (int i = 0; i < n_rows; ++i) {
      Psurf_mat.col(i) = PSurfaceH0(D_.row(i), y_lower, y_upper, sigma, maxP);
    }
  } else {
    double minB = xSeq[0];
    double maxB = xSeq[xSeq.n_elem - 1];
    double dx = (maxB - minB) / resB;
    vec x_lower = seq_C(minB, maxB - dx, dx);
    vec x_upper = x_lower + dx;
    for (int i = 0; i < n_rows; ++i) {
      Psurf_mat.col(i) = PSurfaceHk(D_.row(i), y_lower, y_upper, x_lower, x_upper, sigma, maxP);
    }
  }

  vec out = sum(Psurf_mat, 1);
  return NumericVector(out.begin(),out.end());
}
