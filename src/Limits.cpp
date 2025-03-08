#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
NumericVector findMinMax(const mat& pd, const int& homDim) {
  // Filter rows of dimension homDim
  uvec indices = find(pd.col(0) == homDim); // Get indices of rows where pd.col(0) == homDim

  if (indices.n_elem > 0) {
    // Extract rows with the specified homDim
    vec x = pd.submat(indices, uvec{1});      // Select column 1 elements for matching rows
    vec y = pd.submat(indices, uvec{2});      // Select column 2 elements for matching rows

    // Remove entries with non-finite death times
    uvec finiteIdx = find_finite(y);
    x = x.elem(finiteIdx); // birth
    y = y.elem(finiteIdx); // death
    vec p = y-x; // persistence

    return NumericVector::create(min(x),max(x),max(y),min(p),max(p));
  } else {
    return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL);
  }
}

// [[Rcpp::export]]
NumericVector computeLimits(const field<mat>& Dlist, const int& homDim) {
// Compute extreme values of birth, death and persistence
    int n = Dlist.n_elem;
  NumericVector minB(n), maxB(n), maxD(n), minP(n), maxP(n);

  for (int k = 0; k < n; ++k) {
    NumericVector ans = findMinMax(Dlist(k), homDim);
    minB[k] = ans[0];
    maxB[k] = ans[1];
    maxD[k] = ans[2];
    minP[k] = ans[3];
    maxP[k] = ans[4];
  }

  return NumericVector::create(
    Named("minB") = min(minB),
    Named("maxB") = max(maxB),
    Named("maxD") = max(maxD),
    Named("minP") = min(minP),
    Named("maxP") = max(maxP)
  );
}
