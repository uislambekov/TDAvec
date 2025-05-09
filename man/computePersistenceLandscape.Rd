% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{computePersistenceLandscape}
\alias{computePersistenceLandscape}
\title{Vector Summaries of the Persistence Landscape Functions}

\description{
For a given persistence diagram \eqn{D=\{(b_i,d_i)\}_{i=1}^N} (corresponding to a specified homological dimension), \code{computePersistenceLandscape()} vectorizes the first \eqn{k} persistence landscape functions
\deqn{\lambda_j(t) = j\hbox{max}_{1\leq i \leq N} \Lambda_i(t), \quad j=1,\dots,k,}
where \eqn{j\hbox{max}} returns the \eqn{j}th largest value and
\deqn{\Lambda_i(t) = \left\{
        \begin{array}{ll}
            t-b_i & \quad t\in [b_i,\frac{b_i+d_i}{2}] \\
            d_i-t & \quad t\in (\frac{b_i+d_i}{2},d_i]\\
            0 & \quad \hbox{otherwise}
        \end{array}
    \right.}
based on a scale sequence \code{scaleSeq}. For \emph{generalized} persistence landscape functions, we instead take
\deqn{\Lambda_i(t) = \left\{
        \begin{array}{ll}
            \frac{d_i-b_i}{2K_h(0)}K_h(t-\frac{b_i+d_i}{2}) & \hbox{for } |\frac{t-\frac{b_i+d_i}{2}}{h}| \leq 1 \\
            0 & \hbox{otherwise}
        \end{array}
    \right.}
where \eqn{K_h(\cdot)} is a kernel function with the bandwidth parameter \eqn{h.}
}


\usage{
computePersistenceLandscape(D, homDim, scaleSeq, k = 1, generalized = FALSE, kernel = "triangle", h = NULL)
}
\arguments{
\item{D}{a persistence diagram: a matrix with three columns containing the homological dimension, birth and death values respectively.}

\item{homDim}{the homological dimension (0 for \eqn{H_0}, 1 for \eqn{H_1}, etc.). Rows in \code{D} are filtered based on this value.}

\item{scaleSeq}{a numeric vector of increasing scale values used for vectorization.}

  \item{k}{an integer specifying the number of persistence landscape functions to consider (default is 1).}

  \item{generalized}{a logical value indicating whether to use a generalized persistence landscape or not. Default is FALSE.}

  \item{kernel}{a string specifying the kernel type to use if \code{generalized = TRUE}. Options are "triangle", "epanechnikov", or "tricubic". Default is "triangle".}

  \item{h}{a positive numeric value specifying the bandwidth for the kernel. Must be provided if  \code{generalized = TRUE}.}
}

\details{
The function extracts rows from \code{D} where the first column equals \code{homDim}, and computes values based on the filtered data and \code{scaleSeq}. If \code{D} does not contain any points corresponding to \code{homDim}, a vector of zeros is returned. If \code{generalized = TRUE}, three different kernel functions are currently supported:
  \itemize{
    \item Triangle kernel: \eqn{K_h(t) = \frac{1}{h} \max(0, 1 - \frac{|t|}{h})}
    \item Epanechnikov kernel: \eqn{K_h(t) = \frac{3}{4h} \max(0, 1 - \frac{t^2}{h^2})}
    \item Tricubic kernel: \eqn{K_h(t) = \frac{70}{81h} \max(0, (1 - \frac{|t|^3}{h^3})^3)}
  }
}

\value{
A matrix where the \eqn{j}th column contains the values of the \eqn{j}th order persistence landscape function evaluated at each point of \code{scaleSeq}=\eqn{\{t_1,t_2,\ldots,t_n\}}:

\deqn{
\begin{bmatrix}
\lambda_1(t_1) & \lambda_2(t_1) &  \cdots & \lambda_k(t_1) \\
\lambda_1(t_2) & \lambda_2(t_2) & \cdots & \lambda_k(t_2)\\
\vdots & \vdots& \ddots & \vdots \\
\lambda_1(t_n) & \lambda_2(t_n) & \cdots & \lambda_k(t_n)
\end{bmatrix}
}
}

\author{Umar Islambekov}

\references{
1. Bubenik, P. (2015). Statistical topological data analysis using persistence landscapes. Journal of Machine Learning Research, 16(1), 77-102.

2. Chazal, F., Fasy, B. T., Lecci, F., Rinaldo, A., & Wasserman, L. (2014, June). Stochastic convergence of persistence landscapes and silhouettes. In Proceedings of the thirtieth annual symposium on Computational geometry (pp. 474-483).

3. Berry, E., Chen, Y. C., Cisewski-Kehe, J., & Fasy, B. T. (2020). Functional summaries of persistence diagrams. Journal of Applied and Computational Topology, 4(2):211–262.
}

\examples{
N <- 100 # The number of points to sample

set.seed(123) # Set a random seed for reproducibility

# Sample N points uniformly from the unit circle and add Gaussian noise
theta <- runif(N, min = 0, max = 2 * pi)
X <- cbind(cos(theta), sin(theta)) + rnorm(2 * N, mean = 0, sd = 0.2)

# Compute the persistence diagram using the Rips filtration built on top of X
# The 'threshold' parameter specifies the maximum distance for building simplices
D <- TDAstats::calculate_homology(X, threshold = 2)

scaleSeq = seq(0, 2, length.out = 11) # A sequence of scale values

# Compute vector summaries of the first three persistence landscape functions for homological dimension H_1
computePersistenceLandscape(D, homDim = 1, scaleSeq, k = 3)

# Compute vector summaries of the first three generalized persistence landscape functions (with triangle kernel and bandwidth h=0.2) for homological dimension H_1
computePersistenceLandscape(D, homDim = 1, scaleSeq, generalized = TRUE, k = 3, h = 0.2)

}
