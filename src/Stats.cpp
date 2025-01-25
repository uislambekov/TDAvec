#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Helper function to compute descriptive statistics
vec calcStats(const vec& values) {
  vec P = {0, 0.1, 0.25, 0.5, 0.75, 0.9, 1};
  vec qt = quantile(values, P);
  vec result(9);
  result(0) = mean(values);
  result(1) = stddev(values);
  result(2) = qt(3); // Median
  result(3) = qt(4) - qt(2); // IQR
  result(4) = qt(6) - qt(0); // Range
  result(5) = qt(1); // 10th percentile
  result(6) = qt(2); // 25th percentile
  result(7) = qt(4); // 75th percentile
  result(8) = qt(5); // 90th percentile
  return result;
}

// [[Rcpp::export]]
NumericVector computeStats(const mat& D, const int& homDim) {

  uvec indices = find(D.col(0) == homDim); // Get indices of rows where D.col(0) == homDim

  // If there are no matching rows, stop
  if (indices.n_elem == 0) stop("The diagram has no points corresponding to homological dimension " + std::to_string(homDim));

  // Extract rows with the specified homDim
  vec x = D.submat(indices, uvec{1});      // Select column 1 elements for matching rows
  vec y = D.submat(indices, uvec{2});      // Select column 2 elements for matching rows

  // Total number of bars
  int total_bars = x.n_elem;

  // Remove entries with non-finite death times
  uvec finiteIdx = find_finite(y);
  x = x.elem(finiteIdx);
  y = y.elem(finiteIdx);

  // If x has length zero, stop
  if (x.n_elem == 0) stop("The diagram has no points with finite death value corresponding to homological dimension " + std::to_string(homDim));

  // Compute statistics for births, deaths, midpoints, and lifespans
  vec stats_births = calcStats(x);
  vec stats_deaths = calcStats(y);
  vec stats_midpoints = calcStats((x + y) / 2);
  vec l = y - x; // Lifespans
  vec stats_lifespans = calcStats(l);

  // Entropy calculation
  double L = sum(l);
  double stats_entropy = -sum((l / L) % log2(l / L));

  return NumericVector::create(
    Named("mean_births") = stats_births(0),
    Named("stddev_births") = stats_births(1),
    Named("median_births") = stats_births(2),
    Named("iqr_births") = stats_births(3),
    Named("range_births") = stats_births(4),
    Named("p10_births") = stats_births(5),
    Named("p25_births") = stats_births(6),
    Named("p75_births") = stats_births(7),
    Named("p90_births") = stats_births(8),

    Named("mean_deaths") = stats_deaths(0),
    Named("stddev_deaths") = stats_deaths(1),
    Named("median_deaths") = stats_deaths(2),
    Named("iqr_deaths") = stats_deaths(3),
    Named("range_deaths") = stats_deaths(4),
    Named("p10_deaths") = stats_deaths(5),
    Named("p25_deaths") = stats_deaths(6),
    Named("p75_deaths") = stats_deaths(7),
    Named("p90_deaths") = stats_deaths(8),

    Named("mean_midpoints") = stats_midpoints(0),
    Named("stddev_midpoints") = stats_midpoints(1),
    Named("median_midpoints") = stats_midpoints(2),
    Named("iqr_midpoints") = stats_midpoints(3),
    Named("range_midpoints") = stats_midpoints(4),
    Named("p10_midpoints") = stats_midpoints(5),
    Named("p25_midpoints") = stats_midpoints(6),
    Named("p75_midpoints") = stats_midpoints(7),
    Named("p90_midpoints") = stats_midpoints(8),

    Named("mean_lifespans") = stats_lifespans(0),
    Named("stddev_lifespans") = stats_lifespans(1),
    Named("median_lifespans") = stats_lifespans(2),
    Named("iqr_lifespans") = stats_lifespans(3),
    Named("range_lifespans") = stats_lifespans(4),
    Named("p10_lifespans") = stats_lifespans(5),
    Named("p25_lifespans") = stats_lifespans(6),
    Named("p75_lifespans") = stats_lifespans(7),
    Named("p90_lifespans") = stats_lifespans(8),

    Named("total_bars") = total_bars,
    Named("entropy") = stats_entropy
  );
}
