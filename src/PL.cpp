#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
mat computePL(const mat& D, const int& homDim, const vec& scaleSeq, const int& k=1,
              const bool& generalized=false, const string& kernel = "triangle",
              Nullable<double> h = R_NilValue) {
// Persistence Landscape
  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  mat pl(k,scaleSeq.n_elem);

  // If there are no matching rows, return a zero vector with length scaleSeq.size() - 1
  if (indices.n_elem == 0) return pl.t();

  // If k > number of points in the diagram, then stop
  if (k>indices.n_elem) stop("k must be less than or equal to the number of points in the diagram!");

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  vec Lambda;

  if (generalized){
    double bandwidth;
    if (h.isNotNull()){
      bandwidth = Rcpp::as<double>(h);
      if (bandwidth<0) stop("The bandwidth argument h must be positive!");
    } else stop("Provide the value of the bandwidth argument h");
    vec y_tilde = (y - x) / 2;
    vec x_tilde = (x + y) / 2;
    for (int i=0;i<scaleSeq.n_elem;++i){
      if (kernel == "triangle") {
        Lambda = y_tilde % (1 - (1 / bandwidth) * abs(scaleSeq[i] - x_tilde));
      } else if (kernel == "epanechnikov") {
        Lambda = y_tilde % (1 - (1 / bandwidth) * square(scaleSeq[i] - x_tilde));
      } else if (kernel == "tricubic") {
        Lambda = y_tilde % pow(1 - (1 / bandwidth) * pow(abs(scaleSeq[i] - x_tilde), 3), 3);
      } else {
        stop("Choose between kernel = 'triangle', kernel = 'epanechnikov', or kernel = 'tricubic'");
      }
      Lambda.elem(find(Lambda < 0)).zeros(); // Clamp negative values to 0
      vec sortedLambda = sort(Lambda, "descend");
      uvec kIndices = regspace<uvec>(0, k - 1); // Correct indexing
      pl.col(i) = sortedLambda(kIndices); // Extract the k-th largest elements
    }
  } else{
    for (int i=0;i<scaleSeq.n_elem;++i){
      Lambda = min(scaleSeq[i] - x, y - scaleSeq[i]);
      Lambda.elem(find(Lambda < 0)).zeros(); // Clamp negative values to 0
      vec sortedLambda = sort(Lambda, "descend");
      uvec kIndices = regspace<uvec>(0, k - 1); // Correct indexing
      pl.col(i) = sortedLambda(kIndices); // Extract the k-th largest elements
    }
  }
  return pl.t();
}
