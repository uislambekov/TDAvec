#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericVector computeTropicalCoordinates(const mat& D, const int& homDim, const int&r=1) {

  if (r <= 0) stop("r must be a positive integer!");

  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  // If there are no matching rows, stop
  if (indices.n_elem == 0) stop("The diagram has no points corresponding to homological dimension " + std::to_string(homDim));

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Remove entries with non-finite death times
  uvec finiteIdx = find_finite(y);
  x = x.elem(finiteIdx);
  y = y.elem(finiteIdx);

  // If x has length zero, stop
  if (finiteIdx.n_elem == 0) stop("The diagram has no points with finite death value corresponding to homological dimension " + std::to_string(homDim));

  vec lambda = y - x; //unsorted lifespans
  vec l = sort(y - x,"descend"); // sorted lifespans

  // Compute four most widely-used algebraic functions
  double F1 = l(0);
  double F2;
  double F3;
  double F4;
  if (lambda.n_elem>3){
    F2 = l(0)+l(1);
    F3 = l(0)+l(1)+l(2);
    F4 = l(0)+l(1)+l(2)+l(3);
  } else if (lambda.n_elem==3){
    F2 = l(0)+l(1);
    F3 = l(0)+l(1)+l(2);
    F4 = F3;
  } else if (lambda.n_elem==2){
    F2 = l(0)+l(1);
    F3 = F2;
    F4 = F2;
  }
  double F5 = sum(l);
  vec d = min(r*lambda,x);
  double F6 = sum(d);
  double F7 = sum(max(d+lambda)-(d+lambda));

return NumericVector::create(
  Named("F1") = F1,
  Named("F2") = F2,
  Named("F3") = F3,
  Named("F4") = F4,
  Named("F5") = F5,
  Named("F6") = F6,
  Named("F7") = F7);
}
