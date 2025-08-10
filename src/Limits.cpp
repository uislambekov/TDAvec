#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
vec findMinMax(const mat& pd, const int& homDim) {
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

    return vec({arma::min(x), arma::max(x), arma::max(y), arma::min(p), arma::max(p)});
  } else {
    // Return a vector of NAs if no valid rows are found
    return vec(5, fill::value(arma::datum::nan));
  }
}

// [[Rcpp::export]]
NumericVector computeLimits(const arma::field<arma::mat>& Dlist, const int& homDim) {
// Compute extreme values of birth, death and persistence
    int n = Dlist.n_elem;
  vec minB(n), maxB(n), maxD(n), minP(n), maxP(n);

  for (int k = 0; k < n; ++k) {
    vec ans = findMinMax(Dlist(k), homDim);
    minB[k] = ans[0];
    maxB[k] = ans[1];
    maxD[k] = ans[2];
    minP[k] = ans[3];
    maxP[k] = ans[4];
  }

  uvec finiteIdx = find_finite(minB);

  if (finiteIdx.n_elem==0)
    stop("The diagrams do not contain points corresponding to homological dimension " + std::to_string(homDim));
  else
    return NumericVector::create(
      Named("minB") = min(minB.elem(finiteIdx)),
      Named("maxB") = max(maxB.elem(finiteIdx)),
      Named("maxD") = max(maxD.elem(finiteIdx)),
      Named("minP") = min(minP.elem(finiteIdx)),
      Named("maxP") = max(maxP.elem(finiteIdx))
    );
}
