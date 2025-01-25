#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computeVPB(const mat& D, const int& homDim,
                     const vec& xSeq, const vec& ySeq,
                     const double& tau = 0.3) {
// Vectorized Persistence Block

  // Filter rows matching the given dimension
  uvec indices = find(D.col(0) == homDim);

  // If there are no matching rows, return a zero vector
  if (indices.n_elem == 0) return NumericVector((xSeq.n_elem - 1)*(ySeq.n_elem - 1));

  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Compute lambda and differences
  vec lambda = tau * y;
  vec dy = diff(ySeq);
  int m = dy.n_elem;

  if ((homDim == 0) && (sum(abs(diff(x))) == 0)) {
    vec vpb(m);

    for (int i = 0; i < m; ++i) {
      double c = ySeq[i];
      double d = ySeq[i + 1];

      for (size_t j = 0; j < y.n_elem; ++j) {
        if ((y[j] > c - lambda[j]) && (y[j] < d + lambda[j])) {
          double y_cd = y[j];
          double lambda_cd = lambda[j];
          double yMin = std::max(c, y_cd - lambda_cd);
          double yMax = std::min(d, y_cd + lambda_cd);
          vpb[i] += 0.5 * (yMax * yMax - yMin * yMin) / dy[i];
        }
      }
    }
    return NumericVector(vpb.begin(),vpb.end());
  } else {
    vec dx = diff(xSeq);
    int n = dx.n_elem;
    vec vpb(n * m);

    for (int i = 0; i < m; ++i) {
      double c = ySeq[i];
      double d = ySeq[i + 1];

      uvec yInd = find((y > c - lambda) && (y < d + lambda));
      if (!yInd.is_empty()) {
        vec y_cd = y.elem(yInd);
        vec x_cd = x.elem(yInd);
        vec lambda_cd = lambda.elem(yInd);
        vec yMin_cd = clamp(y_cd-lambda_cd,c,datum::inf);
        vec yMax_cd = clamp(y_cd+lambda_cd,-datum::inf,d);

        for (int j = 0; j < n; ++j) {
          double a = xSeq[j];
          double b = xSeq[j + 1];

          uvec xInd = find((x_cd > a - lambda_cd) && (x_cd < b + lambda_cd));
          if (!xInd.is_empty()) {
            vec x_abcd = x_cd.elem(xInd);
            vec lambda_abcd = lambda_cd.elem(xInd);
            vec xMin = clamp(x_abcd - lambda_abcd,a,datum::inf);
            vec xMax = clamp(x_abcd + lambda_abcd,-datum::inf,b);
            vec yMin = yMin_cd.elem(xInd);
            vec yMax = yMax_cd.elem(xInd);

            vpb[j + i * n] = 0.5 * sum((xMax - xMin) % (yMax - yMin) % (xMin + xMax + yMin + yMax)) / (dx[j] * dy[i]);
          }
        }
      }
    }
    return NumericVector(vpb.begin(),vpb.end());
  }
}
