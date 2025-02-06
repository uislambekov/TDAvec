#include <RcppArmadillo.h>
#include <complex>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// Function S(x,y)
cx_vec S(const vec& x, const vec& y) {
  vec alpha = sqrt(x % x + y % y); // Element-wise norm
  vec factor = (y - x) / (alpha * sqrt(2.0)); // Compute the real factor
  factor.elem(find_nonfinite(alpha)).fill(0); // Set factor to 0 where alpha == 0
  return cx_vec(factor % x, factor % y); // Return complex vector
}

// Function T(x,y)
cx_vec T(const vec& x, const vec& y) {
  vec alpha = sqrt(x % x + y % y); // Element-wise norm
  vec cos_alpha = cos(alpha);     // Element-wise cos(alpha)
  vec sin_alpha = sin(alpha);     // Element-wise sin(alpha)

  vec factor = (y - x) / 2.0; // Compute the real factor
  return cx_vec(factor % (cos_alpha - sin_alpha),
                factor % (cos_alpha + sin_alpha)); // Return complex vector
}

// [[Rcpp::export]]
mat computeComplexPolynomial(const mat& D, const int& homDim, const int& m=1, string polyType = "R") {

    // Get indices of rows where D.col(0) == homDim
  uvec indices = find(D.col(0) == homDim);

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
  if (x.n_elem == 0) stop("The diagram has no points with finite death value corresponding to homological dimension " + std::to_string(homDim));

  // Stop if the number of persistence points is less than m
  if (x.n_elem < m) stop("m must be less than or equal to the number of points in the diagram!");

  // Number of persistence points
  int n = x.n_elem;

  // Initialize a vector to store the coefficients of the polynomial
  cx_vec polynomialCoefficients = {1.0};  // Start with the constant polynomial: P(x) = 1

  // Compute the complex roots depending on polynomial type
  cx_vec roots;
  if (polyType=="R") roots = cx_vec(x, y);
    else if (polyType=="S") roots = S(x,y);
    else if (polyType=="T") roots = T(x,y);
    else stop("Choose between polyType = 'R', polyType = 'S' or polyType = 'T'.");

  // Create the polynomial by multiplying the terms (x - root) for each root
  for (int i = 0; i < n; i++) {
    // Polynomial (x - root_i) * polynomialCoefficients
    cx_vec newPoly = {1.0, -roots(i)};
    polynomialCoefficients = conv(polynomialCoefficients, newPoly);
  }

  // Create the output matrix with real and imaginary parts
  mat output(m, 2);
  output.col(0) = real(polynomialCoefficients.subvec(1, m));  // Real part
  output.col(1) = imag(polynomialCoefficients.subvec(1, m));  // Imaginary part

  return output;
}
