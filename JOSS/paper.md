---
title: 'TDAvec: Computing  Vector Summaries of Persistence Diagrams for Topological
  Data Analysis in R and Python'
authors:
- name: Umar Islambekov
  orcid: "0000-0002-5103-3238"
  affiliation: 1
  equal-contrib: yes
- name: Aleksei Luchinsky
  email: aluchi@bgsu.edu
  affiliation: 1
  corresponding: yes
date: "2025-04-14"
output: pdf_document
bibliography: paper.bib
tags:
- Persistent Homology
- Topological Data Analysis
- Persistence Diagram
- Vectorization of Persistence Diagrams
affiliations:
- index: 1
  name: Bowling Green State University, USA
---


# Summary

The theory of \emph{persistent homology} is one of the popular tools in \emph{topological data analysis} (TDA) to analyze data with underlying shape structure [@Carlsson:2009; @edelsbrunner2010computational; @chazal2021introduction]. In this context, a single data observation could be a collection of points lying in a metric space, an image, a graph or a time series. The basic idea behind persistent homology is to build a nested sequence (or \emph{filtration}) of \emph{simplicial complexes} (indexed by a scale parameter) on top of data points and keep a record of the appearance and disappearance of various topological features at different scale values. Here, these topological features are "holes" of different dimensions -- connected components, loops, voids, and their higher-dimensional versions whose emergence and subsequent disappearance are tracked using a concept of homology from algebraic topology. From a geometric point of view, simplicial complexes consist of vertices, edges, triangles, tetrahedra etc., glued together and serve as a means for recovering (at least partially) the underlying shape information which is lost during sampling [@nanda2013simplicial].

A topological descriptor outputted by the persistent homology encoding the shape of data is called a \emph{persistence diagram} (PD). Mathematically, a $k$-dimensional PD is a multi-set of points $D=\{(b_i,d_i)\}_{i=1}^N$, where each point $(b_i,d_i)$ corresponds to a topological feature of homological dimension $k$ (0 if a connected component, 1 if a loop, 2 if a void, etc) with the $x$-coordinate representing the scale at which this feature is born (or created), and the $y$-coordinate representing the scale at which it dies (or disappears). In practice, one is usually interested in applying a machine learning method to PDs to make further inferences from data. However, the fact that PDs do not form a Hilbert space, which is a feature (or an input) space for a wide class of machine learning methods, limits their direct of use in applications. To overcome this challenge, kernel methods and vectorization techniques are commonly used [@chung2022persistence]. The kernel approach involves defining a notion of similarity between pairs of PDs, whereas the vectorization methods aim to transform PDs into finite-dimensional feature vectors that can be used as input for many standard machine learning models. In recent years, the kernel and vectorization approaches have proven successful and gained  prominence in the applied TDA literature (see @hensel2021survey for a survey of applications of TDA in machine learning).   

### R tools for TDA
The computational tools for TDA in the \texttt{R} environment are provided through various packages\footnote{In this overview, we only focus on \texttt{R} packages for TDA that are available on the CRAN repository.} such as \texttt{TDA} [@TDA], \texttt{TDAstats} [@wadhwa2018tdastats], \texttt{kernelTDA} [@kernelTDA], \texttt{TDAmapper} [@TDAmapper], \texttt{TDAkit} [@TDAkit], \texttt{tdaunif} [@tdaunif], \texttt{TDApplied} [@TDApplied] and \texttt{ripserr} [@ripserr]. 

The \texttt{TDA} package is the largest \texttt{R} package for TDA. The \texttt{TDA} package offers tools to compute PDs for commonly used types of filtrations such as \emph{Vietoris-Rips}, \emph{Alpha} and \emph{Alpha shape}. It also allows to construct more general sublevel set filtrations and compute the corresponding PDs. Moreover, the \texttt{TDA} package provides implementations to plot PDs and compute \emph{bottleneck and Wasserstein} distances between them. \texttt{TDAstats} offers a variety of tools for conducting statistical inference (such as hypothesis testing) on PDs. Compared to the \texttt{TDA} package, it computes PDs much faster for Vietoris-Rips filtrations based on the Ripser C++ library [@Bauer2021Ripser] and offers more aesthetic visualization of the diagrams using the \texttt{ggplot2} package [@ggplot2]. The \texttt{kernelTDA} package contains implementations of popular kernel-based methods for TDA such as \emph{geodesic Gaussian kernel}, \emph{geodesic Laplacian kernel}, \emph{persistence Fisher kernel} and \emph{persistence sliced Wasserstein kernel}. For computing the Wasserstein distance between a pair of PDs, unlike the \texttt{TDA} package, it uses an iterative procedure to reasonably approximate the exact distance which that leads to a considerable reduction in run-time cost. The \texttt{ripserr} package allows a fast computation of PDs for filtrations on Vietoris-Rips and cubical complexes using the Ripser C++ library. The \texttt{TDApplied} and \texttt{TDAkit} packages provides various tools to integrate topological features (PDs or their vector summaries) into machine and statistical learning settings. The \texttt{tdaunif} is a useful package if one needs to sample points from various manifolds such as a klein-bottle, an ellipse or a torus. \texttt{TDAmapper} offers tools to visualize high-dimensional data by constructing the so-called Mapper graphs that preserve its topological structure.

### Python tools for TDA

Several \texttt{Python} libraries are available for TDA, such as \texttt{Giotto-tda} [@tauzin2021giotto], \texttt{Ripser} [@christopher2018lean], \texttt{Gudhi} [@rouvreau2020gudhi], \texttt{Scikit-tda}, \texttt{Dionysus 2} [@Dionysus2], \texttt{Persim} [@Persim] and \texttt{KeplerMapper} [@van2019kepler]. \texttt{Giotto-tda} is a powerful library that integrates with the popular machine learning library \texttt{scikit-learn}, offering tools for persistent homology and visualizations of persistence diagrams. \texttt{Ripser} focuses on fast computation of Vietoris-Rips complexes, especially for large datasets. \texttt{Gudhi} provides a wide range of topological tools for simplicial complexes, persistent homology, and topological signatures. \texttt{Scikit-tda} is another package that integrates with \texttt{scikit-learn}, simplifying the application of TDA to typical machine learning tasks. \texttt{Dionysus 2} offers fast computation of persistent homology and cohomology, with an emphasis on flexibility and efficiency. \texttt{Persim} focuses on tools for working with PDs. It contains implementations of commonly used vectorization and kernel methods for PDs. \texttt{KeplerMapper} implements the TDA Mapper algorithm to visualize high-dimensional data. For a more comprehensive list of \texttt{Python} libraries for TDA and their functionality, we refer the readers to [@awesome-tda:2024].

# Statement of need

The problem of transforming PDs into finite dimensional vectors for machine learning purposes has attracted considerable attention in the TDA research community over the past decade. While early vectors summaries of PDs such as \emph{persistence landscape} [@bubenik2015statistical], \emph{persistence silhouette} [@chazal2014stochastic], \emph{Betti curve}\footnote{also called \emph{Betti function}} [@chazal2021introduction] and \emph{persistence image} [@adams2017persistence] have been implemented in both \texttt{Python} and \texttt{R} packages, there is no single package that systematically gathers them - the early ones as well as those introduced in recent years - in one place under unified syntax. Moreover, in \texttt{R}, all the code behind the existing vector implementations is written using standard \texttt{R} functions which may prove slow and inefficient for large-scale computations. The \texttt{TDAvec} \texttt{R} package and its \texttt{Python} version aim to fill in some of these gaps. Its contributions can be summarized in the following three areas: 

1. It expands the list of implemented vector summaries for PDs by providing vectorizations of 13 commonly used methods in TDA. These methods are grouped into three broad categories:

  - Functional vector summaries - based on summary functions:
    - Betti curve [@chazal2021introduction]
    - Euler characteristic curve [@richardson2014efficient]
    - Normalized life curve [@chung2022persistence]
    - Persistence block  [@chan2022computationally]
    - Persistence surface [@adams2017persistence]
    - Persistence landscape function [@bubenik2015statistical]
    - Persistence silhouette function [@chazal2014stochastic]
    - Persistent entropy summary function [@atienza2020stability]
    - Template function [@perea2023approximating]

  - Algebraic vector summaries - based on polynomial maps:
    - Algebraic functions [@Algebraic_functions]
    - Complex polynomial coefficients [@ferri1999representing; @di2015comparing]
    - Tropical coordinate function [@kalivsnik2019tropical]
     
  - Statistical vector summaries - based on descriptive statistics:
    - Basic descriptive statistics [@ali2023survey]

2. A univariate summary function $f$ of a PD is commonly vectorized by evaluating it at a sequence of points on a one-dimensional grid, then organizing the resulting values into a vector:
\begin{equation}\label{stand_vec} (f(t_1),f(t_2),\ldots,f(t_n))\in\mathbb{R}^n, \end{equation}
where $t_1, t_2, \ldots, t_n$ form an increasing sequence of scale values. For instance, the \texttt{landscape()} and \texttt{silhouette()} functions in the \texttt{TDA} package produce such vector summaries for persistence landscapes and silhouettes, respectively.
In addition to this standard approach, the \texttt{TDAvec} package introduces an alternative vectorization scheme that captures the average behavior of $f$ between consecutive scale values $t_i$ and $t_{i+1}$ through integration:
\begin{equation} \Big(\frac{1}{\Delta t_1}\int_{t_1}^{t_2}f(t),dt,\frac{1}{\Delta t_2}\int_{t_2}^{t_3}f(t),dt,\ldots,\frac{1}{\Delta t_{n-1}}\int_{t_{n-1}}^{t_n}f(t),dt\Big)\in\mathbb{R}^{n-1}, \end{equation}
where $\Delta t_i = t_{i+1} - t_i$. Unlike the method in (\ref{stand_vec}), this approach retains information about the behavior of $f$ between neighboring scale points. It is applicable to any univariate summary function that is integrable in closed form, such as the persistence silhouette, persistent entropy summary function, Euler characteristic curve, normalized life curve, and Betti function. Users have the flexibility to choose between the two vectorization methods based on their application needs.

3. To achieve higher computational efficiency, all code behind the vector summaries of \texttt{TDAvec} is written in C++ using the `Rcpp` and `RcppArmadillo` packages. 

The \texttt{TDAvec} \texttt{R} package and a vignette showing its basic usage with examples are available on the CRAN repository\footnote{https://cran.r-project.org/web/packages/TDAvec/index.html}. For \texttt{Python} examples, we refer the readers to https://github.com/ALuchinsky/tdavec/.

# References
