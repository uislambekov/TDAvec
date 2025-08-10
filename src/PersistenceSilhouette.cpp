#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector computePersistenceSilhouette(const arma::mat& D, const int& homDim, const arma::vec& scaleSeq,
                        const double& p=1.0, const std::string& evaluate = "intervals") {

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

  // Compute weights based on powered differences
  vec pp = pow(y - x, p);
  vec w = pp / sum(pp);

  // Prepare output vector
  vec phi;

  if (evaluate=="intervals"){
    // Initialize output vector phi
    int l = scaleSeq.n_elem - 1;
    phi.set_size(l);

    // Compute the phi values
    for (int k = 0; k < l; ++k) {
      vec alpha1 = clamp(x,scaleSeq[k],datum::inf); //pmax(x,scaleSeq[k])
      vec alpha2 = clamp((x + y) / 2.0,scaleSeq[k],datum::inf);
      vec beta1 = clamp((x + y) / 2.0,-datum::inf,scaleSeq[k+1]);
      vec beta2 = clamp(y,-datum::inf,scaleSeq[k+1]); //pmin(y,scaleSeq[k+1])

      vec b1 = clamp(beta1 - alpha1,0,datum::inf) % ((beta1 + alpha1) / 2.0 - x);
      vec b2 = clamp(beta2 - alpha2,0,datum::inf) % (y - (beta2 + alpha2) / 2.0);

      phi[k] = sum(w % (b1 + b2)) / (scaleSeq[k+1] - scaleSeq[k]);
    }
  } else if (evaluate=="points") {
    // Initialize output vector phi
    int l = scaleSeq.n_elem;
    phi.set_size(l);

    // Compute the phi values
    for (int k = 0; k < l; ++k) {
      vec b1 = clamp(scaleSeq[k]-x,0,datum::inf);
      vec b2 = clamp(y-scaleSeq[k],0,datum::inf);
      phi[k] = sum(w % min(b1,b2)); // min(b1,b2) is pmin(b1,b2) in R
    }
  }

  return NumericVector(phi.begin(),phi.end());
}
