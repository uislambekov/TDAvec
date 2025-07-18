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
date: "2025-06-16"
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

A topological descriptor outputted by the persistent homology encoding the shape of data is called a \emph{persistence diagram} (PD). Mathematically, a $k$-dimensional PD is a multi-set of points $D=\{(b_i,d_i)\}_{i=1}^N$, where each point $(b_i,d_i)$ corresponds to a topological feature of homological dimension $k$ (0 if a connected component, 1 if a loop, 2 if a void, etc) with the $x$-coordinate representing the scale at which this feature is born (or created), and the $y$-coordinate representing the scale at which it dies (or disappears). In practice, one is usually interested in applying a machine learning method to PDs to make further inferences from data. However, the fact that PDs do not form a Hilbert space, which is a feature (or an input) space for a wide class of machine learning methods, limits their direct of use in applications. To overcome this challenge, kernel methods and vectorization techniques are commonly used [@chung2022persistence]. The kernel approach involves defining a notion of similarity between pairs of PDs, whereas the vectorization methods aim to transform PDs into finite-dimensional feature vectors that can be used as input for many standard machine learning models.   

# Statement of need

The problem of transforming PDs into finite dimensional vectors for machine learning purposes has attracted considerable attention in the TDA research community over the past decade. Early vector summaries of PDs - such as the persistence landscape [@bubenik2015statistical], persistence silhouette [@chazal2014stochastic], Betti curve\footnote{Also known as the Betti function} [@chazal2021introduction], and persistence image [@adams2017persistence] - have been implemented in both \texttt{Python} and \texttt{R} packages. However, there remains a need for a unified package that brings these methods (both the classical approaches and those developed more recently) together using consistent syntax and efficient implementation. The \texttt{TDAvec} package, available in both \texttt{R} and \texttt{Python}, is designed to meet this need. Its contributions can be summarized in the following three areas: 

1. It extends the list of implemented vector summaries for PDs by incorporating 13 vectorization methods used in TDA. These methods can be grouped into three broad categories:

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

The \texttt{TDAvec} \texttt{R} package and a vignette showing its basic usage with examples are available on the CRAN repository\footnote{https://cran.r-project.org/web/packages/TDAvec/index.html}. For \texttt{Python} examples, we refer the readers to [this page](https://github.com/ALuchinsky/tdavec/).

# References
