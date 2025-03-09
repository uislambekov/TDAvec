#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computeAlgebraicFunctions(const mat& D, const int& homDim) {

  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  // If there are no matching rows, return zeros
  if (indices.n_elem == 0) return NumericVector::create(
    Named("f1") = 0,
    Named("f2") = 0,
    Named("f3") = 0,
    Named("f4") = 0);

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Remove entries with non-finite death times
  uvec finiteIdx = find_finite(y);
  x = x.elem(finiteIdx);
  y = y.elem(finiteIdx);

  // If x has length zero, return zeros
  if (x.n_elem == 0) return NumericVector::create(
    Named("f1") = 0,
    Named("f2") = 0,
    Named("f3") = 0,
    Named("f4") = 0);

  double yMax = max(y);
  vec l = y - x;
  vec lLeft = yMax-y;
  vec x2 = pow(x,2);
  vec lLeft2 = pow(lLeft,2);
  vec l4 = pow(l,4);

  // Compute four most widely-used algebraic functions
  double f1 = sum(x % l);
  double f2 = sum(lLeft % l);
  double f3 = sum(x2 % l4);
  double f4 = sum(lLeft2 % l4);

  return NumericVector::create(
  Named("f1") = f1,
  Named("f2") = f2,
  Named("f3") = f3,
  Named("f4") = f4);
}
