% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{computePersistenceImage}
\alias{computePersistenceImage}
\title{A Vector Summary of the Persistence Surface}

\description{
For a given persistence diagram \eqn{D=\{(b_i,p_i)\}_{i=1}^N} (corresponding to a specified homological dimension), \code{computePersistenceImage()} computes the persistence image - a vector summary of the persistence surface: \deqn{\rho(x,y)=\sum_{i=1}^N f(b_i,p_i)\phi_{(b_i,p_i)}(x,y),} where \eqn{\phi_{(b_i,p_i)}(x,y)} is
the Gaussian distribution with mean \eqn{(b_i,p_i)} and
covariance matrix \eqn{\sigma^2 I_{2\times 2}} and
\deqn{f(b,p) = w(p)=\left\{
        \begin{array}{ll}
            0 & \quad p\leq 0 \\
            p/p_{max} & \quad 0<p<p_{max}\\
            1& \quad p\geq p_{max}
        \end{array}
    \right.}
is the weighting function with \eqn{p_{max}} being the maximum persistence value among all persistence diagrams considered in the experiment. Points in \eqn{D} with infinite persistence values are ignored.
}

\usage{
computePersistenceImage(D, homDim, xSeq, ySeq, sigma)
}

\arguments{
\item{D}{a persistence diagram: a matrix with three columns containing the homological dimension, \emph{birth} and \emph{persistence} values respectively.}

\item{homDim}{the homological dimension (0 for \eqn{H_0}, 1 for \eqn{H_1}, etc.). Rows in \code{D} are filtered based on this value.}

\item{xSeq}{a numeric vector of increasing x (birth) values used for vectorization.}

\item{ySeq}{a numeric vector of increasing y (persistence) values used for vectorization.}

\item{sigma}{a standard deviation (\eqn{\sigma}) of the Gaussian.}
}

\details{
The function extracts rows from \code{D} where the first column equals \code{homDim}, and computes values based on the filtered data, \code{xSeq} and \code{ySeq}. If \code{D} does not contain any points corresponding to \code{homDim}, a vector of zeros is returned.
}

\value{
A numeric vector whose elements are the average values of the persistence surface computed over each cell of the two-dimensional grid constructred from \code{xSeq}=\eqn{\{x_1,x_2,\ldots,x_n\}} and \code{ySeq}=\eqn{\{y_1,y_2,\ldots,y_m\}}:

\deqn{\Big(\frac{1}{\Delta x_1\Delta y_1}\int_{[x_1,x_2]\times [y_1,y_2]}\rho(x,y)dA,\ldots,\frac{1}{\Delta x_{n-1}\Delta y_{m-1}}\int_{[x_{n-1},x_n]\times [y_{m-1},y_m]}\rho(x,y)dA\Big)\in\mathbb{R}^{d},}
where \eqn{d=(n-1)(m-1)}, \eqn{dA=dxdy}, \eqn{\Delta x_k=x_{k+1}-x_k} and \eqn{\Delta y_j=y_{j+1}-y_j.}

If \code{homDim=0} and all the birth values are equal (e.g., zero), univariate Gaussians are used instead for vectorization:

\deqn{\Big(\frac{1}{\Delta y_1}\int_{y_1}^{y_2}\rho(y)dy,\ldots,\frac{1}{\Delta y_{m-1}}\int_{y_{m-1}}^{y_m}\rho(y)dy\Big)\in\mathbb{R}^{m-1},}
where \eqn{\rho(y)=\sum_{i=1}^N f(p_i)\phi_{p_i}(y)} and \eqn{\Delta y_j=y_{j+1}-y_j.}

}

\author{Umar Islambekov}

\references{
1. Adams, H., Emerson, T., Kirby, M., Neville, R., Peterson, C., Shipman, P., ... & Ziegelmeier, L. (2017). Persistence images: A stable vector representation of persistent homology. Journal of Machine Learning Research, 18.
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

# Switch from the birth-death to the birth-persistence coordinates
D[,3] <- D[,3] - D[,2]
colnames(D)[3] <- "Persistence"

resB <- 5 # Resolution (or grid size) along the birth axis
resP <- 5 # Resolution (or grid size) along the persistence axis

# Compute persistence image for homological dimension H_0
minPH0 <- min(D[D[,1]==0,3]); maxPH0 <- max(D[D[,1]==0,3])
ySeqH0 <- seq(minPH0, maxPH0, length.out = resP + 1)
sigma <- 0.5 * (maxPH0 - minPH0) / resP
computePersistenceImage(D, homDim = 0, xSeq = NA, ySeq = ySeqH0, sigma = sigma)

# Compute persistence image for homological dimension H_1
minBH1 <- min(D[D[,1]==1,2]); maxBH1 <- max(D[D[,1]==1,2])
minPH1 <- min(D[D[,1]==1,3]); maxPH1 <- max(D[D[,1]==1,3])
xSeqH1 <- seq(minBH1, maxBH1, length.out = resB + 1)
ySeqH1 <- seq(minPH1, maxPH1, length.out = resP + 1)
sigma <- 0.5 * (maxPH1 - minPH1) / resP
computePersistenceImage(D, homDim = 1, xSeq = xSeqH1, ySeq = ySeqH1, sigma = sigma)

}
