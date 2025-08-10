#include <RcppArmadillo.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
double tent_function_2D(const vec& x, const vec& y, const double &a, const double &b, const double& delta) {
  // Compute the max of |x-a| and |y-b|
  vec max_dist = max(abs(x - a),abs(y - b));
  // Compute g_{(a,b),delta}(x)
  double result = sum(clamp(1.0 - (1.0 / delta) * max_dist,0, datum::inf));
  return result;
}

double tent_function_1D(const vec& y, const double& b, const double& delta) {
  // Compute g_{a,delta}(x)
  double result = sum(clamp(1.0 - (1.0 / delta) * abs(y - b),0, datum::inf));
  return result;
}

// [[Rcpp::export]]
NumericVector computeTemplateFunction(const arma::mat& D, const int& homDim, const double& delta, const int& d, const double& epsilon) {
  // delta = increment
  // d = number of bins in each axis
  // epsilon = vertical shift up

  if (delta<0) stop("The argument 'delta' must be positive!");
  if (epsilon<0) stop("The argument 'epsilon' must be positive!");

  // Get indices of rows where D.col(0) == homDim
  uvec indices = find(D.col(0) == homDim);

  // If there are no matching rows, return a zero vector
  if (indices.n_elem == 0) return NumericVector((d+1)*d);

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Remove entries with non-finite death times
  uvec finiteIdx = find_finite(y);
  x = x.elem(finiteIdx);
  y = y.elem(finiteIdx);

  if (any(x < 0)) stop("The birth values must all be positive!");
  if (x.n_elem == 0) return NumericVector((d+1)*d);

  // Compute lifespans
  vec l = y-x;

  // Initialize the output vector
  vec tf;
  double sumX = sum(abs(diff(x)));

  if ((homDim == 0) && (sumX == 0)) {
    tf.set_size(d);
    // Compute center_l
    vec center_l = linspace(delta, d*delta, d)+epsilon;
    for (int j = 0; j < center_l.n_elem; ++j) {
      // Sum the values and store in tf
      tf(j) = tent_function_1D(l,center_l(j),delta);
    }
  } else {
    tf.set_size((d+1)*d);
    // Compute center_x and center_l
    vec center_x = linspace(0, d*delta, d+1);
    vec center_l = linspace(delta, d*delta, d)+epsilon;
    // Index to iterate through the grid
    int idx = 0;
    for (int i = 0; i < center_x.n_elem; ++i) {
      for (int j = 0; j < center_l.n_elem; ++j) {
        // Current grid point (a, b)
        double a = center_x(i);
        double b = center_l(j);
        // Sum the values and store in tf
        tf(idx) = tent_function_2D(x, l, a, b, delta);
        ++idx;
      }
    }
  }
  return NumericVector(tf.begin(), tf.end());  // Convert Armadillo vector to NumericVector
}
