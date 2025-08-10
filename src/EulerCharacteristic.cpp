#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
vec computeBetti(const mat& D, const int& homDim, const vec& scaleSeq, string evaluate = "intervals") {
  // Get indices of rows where D.col(0) == homDim
  uvec indices = find(D.col(0) == homDim);
  // If there are no matching rows, return a zero vector with length scaleSeq.size() - 1
  if (indices.n_elem == 0) {
    if (evaluate=="intervals") {
      return NumericVector(scaleSeq.n_elem - 1);
    } else if (evaluate=="points") {
      return NumericVector(scaleSeq.n_elem);
    }
  }
  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  vec betti;

  if (evaluate=="intervals"){
    // Initialize output vector betti
    int l = scaleSeq.n_elem - 1;
    betti.set_size(l);

    // Compute the betti values
    for (int k = 0; k < l; ++k) {
      // Element-wise operations using pmin and pmax logic
      vec b = clamp(y, -datum::inf, scaleSeq[k+1]) - clamp(x, scaleSeq[k],datum::inf);
      betti[k] = sum(clamp(b, 0, datum::inf)) / (scaleSeq[k+1] - scaleSeq[k]);
    }
  } else if (evaluate=="points"){
    // Initialize output vector betti
    int l = scaleSeq.n_elem;
    betti.set_size(l);

    //Compute the betti values
    for (int k = 0; k < l; ++k) {
      betti[k] = sum((scaleSeq[k] >= x) && (scaleSeq[k] < y));
    }
  }
  return betti;
}

// [[Rcpp::export]]
NumericVector computeEulerCharacteristic(const arma::mat& D, const arma::vec& scaleSeq, const int& maxhomDim = -1, const std::string& evaluate = "intervals") {

  int effectiveMaxHomDim;

  if (maxhomDim == -1) {
    effectiveMaxHomDim = max(D.col(0)); // Max of homological dimensions in column 1
  } else {
    effectiveMaxHomDim = maxhomDim;
  }

  // Initialize an Armadillo matrix for ecc
  mat ecc;

  if (evaluate=="intervals") {
    ecc.set_size(scaleSeq.n_elem - 1, effectiveMaxHomDim + 1);
  } else if (evaluate=="points") {
    ecc.set_size(scaleSeq.n_elem, effectiveMaxHomDim + 1);
  } else stop("Choose between evaluate = 'intervals' or evaluate = 'points'");


  // Loop through each dimension
  for (int d = 0; d <= effectiveMaxHomDim; ++d) {
    vec betti = computeBetti(D,d,scaleSeq,evaluate);
    // Multiply scalar sign with the result of computeBetti
    ecc.col(d) = pow(-1,d) * betti;  // Assign to the column
  }

  // Compute and return the row sums
  vec ans = sum(ecc, 1); // Calculate the row sums
  return NumericVector(ans.begin(),ans.end());
}
