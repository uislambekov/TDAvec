#include <RcppArmadillo.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computePES(const mat& D, const int& homDim, const vec& scaleSeq, string evaluate = "intervals") {
// PES = Persistence Entropy Summary

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

  // Calculate lL (either difference or normalized difference)
  vec lL;
  if (x.n_elem == 1) {
    lL = y - x;  // Just the difference if there's only one row
  } else {
    lL = (y - x) / sum(y - x);  // Normalized difference
  }

  // Calculate entropy values
  vec entr = -lL % log2(lL);  // Element-wise operations

  // Output vector
  vec pes;

  if (evaluate == "intervals") {
    int l = scaleSeq.n_elem - 1;
    pes.set_size(l);

    // Compute the pes values for intervals
    for (int k = 0; k < l; ++k) {
      // Element-wise operations using pmin and pmax logic
      vec b = clamp(y, -datum::inf, scaleSeq[k+1]) - clamp(x, scaleSeq[k], datum::inf);
      pes[k] = sum(entr % clamp(b, 0, datum::inf)) / (scaleSeq[k+1] - scaleSeq[k]);
    }

  } else if (evaluate == "points") {
    int l = scaleSeq.n_elem;
    pes.set_size(l);

    // Compute the pes values for points
    for (int k = 0; k < l; ++k) {
      // Element-wise logical operation on x and y with respect to scaleSeq[k]
      pes[k] = sum(entr % ((scaleSeq[k] >= x) && (scaleSeq[k] < y)));
    }
  } else {
    stop("Choose between evaluate = 'intervals' or evaluate = 'points'");
  }

  return NumericVector(pes.begin(), pes.end());  // Convert Armadillo vector to NumericVector
}
