#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computeBetti(const mat& D, const int& homDim, const vec& scaleSeq, const string& evaluate = "intervals") {

  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  // If there are no matching rows, return a zero vector with length scaleSeq.size() - 1
  if (indices.n_elem == 0) {
    if (evaluate=="intervals") {
      return NumericVector(scaleSeq.n_elem - 1);
    } else if (evaluate=="points") {
        return NumericVector(scaleSeq.n_elem);
    } else stop("Choose between evaluate = 'intervals' or evaluate = 'points'");
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

  } else stop("Choose between evaluate = 'intervals' or evaluate = 'points'");

    return NumericVector(betti.begin(),betti.end());
}
