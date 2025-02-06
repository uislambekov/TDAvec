#include <RcppArmadillo.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computeNormalizedLife(const arma::mat& D, const int& homDim, const arma::vec& scaleSeq, std::string evaluate = "intervals") {

  // Get indices of rows where D.col(0) == homDim
  uvec indices = find(D.col(0) == homDim);

  // If there are no matching rows, return a zero vector
  if (indices.n_elem == 0) {
    if (evaluate == "intervals") {
      return NumericVector(scaleSeq.n_elem - 1);
    } else if (evaluate == "points") {
      return NumericVector(scaleSeq.n_elem);
    } else {
      stop("Choose between evaluate = 'intervals' or evaluate = 'points'");
    }
  }

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Remove entries with non-finite death times
  uvec finiteIdx = find_finite(y);
  x = x.elem(finiteIdx);
  y = y.elem(finiteIdx);

  if (x.n_elem == 0) {
    if (evaluate == "intervals") {
      return NumericVector(scaleSeq.n_elem - 1);
    } else if (evaluate == "points") {
      return NumericVector(scaleSeq.n_elem);
    }
  }

  // Calculate lL (normalized difference)
  vec lL = (y - x) / sum(y - x);

  // Output vector
  vec nl;

  if (evaluate == "intervals") {
    int l = scaleSeq.n_elem - 1;
    nl.set_size(l);

    // Compute the nl values for intervals
    for (int k = 0; k < l; ++k) {
      // Element-wise operations using pmin and pmax logic
      vec b = clamp(y, -datum::inf, scaleSeq[k+1]) - clamp(x, scaleSeq[k], datum::inf);
      nl[k] = sum(lL % clamp(b, 0, datum::inf)) / (scaleSeq[k+1] - scaleSeq[k]);
    }

  } else if (evaluate == "points") {
    int l = scaleSeq.n_elem;
    nl.set_size(l);

    // Compute the nl values for points
    for (int k = 0; k < l; ++k) {
      // Element-wise logical operation on x and y with respect to scaleSeq[k]
      nl[k] = sum(lL % ((scaleSeq[k] >= x) && (scaleSeq[k] < y)));
    }
  } else {
    stop("Choose between evaluate = 'intervals' or evaluate = 'points'");
  }

  return NumericVector(nl.begin(), nl.end());  // Convert Armadillo vector to NumericVector
}
