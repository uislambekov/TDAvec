---
title: 'TDAvec: An R Package for Computing  Vector Summaries of Persistence Diagrams for Topological Data Analysis'
tags:
- R
- Persistent Homology
- Topological Data Analysis
- Persistence Diagram
- Vectorization of Persistence Diagram
date: "8 August 2024"
output: pdf_document
authors:
- name: Umar Islambekov
  orcid: "0000-0002-5103-3238"
  equal-contrib: yes
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Bowling Green State University, USA
  index: 1
---

# Summary

The theory of \emph{persistent homology} is one of the popular tools in \emph{topological data analysis} (TDA) to analyze data with underlying shape structure [@Carlsson:2009; @edelsbrunner2010computational; @chazal2021introduction]. In this context, a single data observation could be a collection of points lying in a metric space, an image, a graph or a time series. The basic idea behind persistent homology is to build a nested sequence (or \emph{filtration}) of \emph{simplicial complexes} (indexed by a scale parameter) on top of data points and keep a record of the appearance and disappearance of various topological features at different scale values. Here, these topological features are "holes" of different dimensions -- connected components, loops, voids, and their higher-dimensional versions whose emergence and subsequent disappearance are tracked using a concept of homology from algebraic topology. From a geometric point of view, simplicial complexes consist of vertices, edges, triangles, tetrahedra etc., glued together and serve as a means for recovering (at least partially) the underlying shape information which is lost during sampling [@nanda2013simplicial].

A topological descriptor outputted by the persistent homology encoding the shape of data is called a \emph{persistence diagram} (PD). Mathematically, a $k$-dimensional PD is a multi-set of points $D=\{(b_i,d_i)\}_{i=1}^N$, where each point $(b_i,d_i)$ corresponds to a topological feature of homological dimension $k$ (0 if a connected component, 1 if a loop, 2 if a void, etc) with the $x$-coordinate representing the scale at which this feature is born (or created), and the $y$-coordinate representing the scale at which it dies (or disappears). In practice, one is usually interested in applying a machine learning method to PDs to make further inferences from data. However, the fact that PDs do not form a Hilbert space, which is a feature (or an input) space for a wide class of machine learning methods, limits their direct of use in applications. To overcome this challenge, kernel methods and vectorization techniques are commonly used [@chung2022persistence]. The kernel approach involves defining a notion of similarity between pairs of PDs, whereas the vectorization methods aim to transform PDs into finite-dimensional feature vectors that can be used as input for many standard machine learning models. Such vector summaries of PDs are computed in two steps: first one constructs a respective summary function from a given PD and then vectorizes it using either one or two dimensional grid of scale values. In recent years, the kernel and vectorization approaches have proven successful and gained  prominence in the applied TDA literature (see @hensel2021survey for a survey of applications of TDA in machine learning).   

The computational tools for TDA in the \texttt{R} environment are provided through various packages such as \texttt{TDA} [@TDA], \texttt{TDAstats} [@wadhwa2018tdastats], \texttt{kernelTDA} [@kernelTDA], \texttt{TDAmapper} [@TDAmapper], \texttt{TDAkit} [@TDAkit], \texttt{tdaunif} [@tdaunif], \texttt{TDApplied} [@TDApplied] and \texttt{ripserr} [@ripserr] (see Table \ref{tab:TDA_packages} for an overview of these packages in terms of their scope and areas of focus). 

\begin{table}[htbp]
  \centering
  \begin{tabular}{|l|c|c|c|c|c|c|c|c|}
    \hline
    & \rotatebox{90}{\tt TDA} & \rotatebox{90}{\tt TDApplied} & \rotatebox{90}{\tt TDAstats} & \rotatebox{90}{\tt TDAkit} & \rotatebox{90}{\tt kernelTDA} &  \rotatebox{90}{\tt tdaunif} &\rotatebox{90}{\tt ripserr} & \rotatebox{90}{\tt TDAmapper} \\
    \hline    
    Sampling methods & \checkmark & & &\checkmark & & \checkmark& &\\
    \hline
    Density estimation & \checkmark & & & & & & &\\
    \hline
    Alpha filtration & \checkmark & & & & & & &\\
    \hline
    Alpha shape filtration &\checkmark & & & & & & &\\
    \hline
    Vietoris-Rips filtration &\checkmark & &\checkmark &\checkmark & & & \checkmark&\\
    \hline
    User-defined filtration &\checkmark & & & & & & &\\
    \hline
    Cubical complex & & & & & & & \checkmark &\\
    \hline
    Wasserstein distance &\checkmark &\checkmark &\checkmark & & \checkmark & & &\\    
    \hline
    Plotting persistence diagrams &\checkmark &\checkmark &\checkmark &\checkmark & & & &\\    
    \hline
    Statistical methods &\checkmark &\checkmark &\checkmark &\checkmark & & & &\\    
    \hline
    Vectorization methods &\checkmark & & &\checkmark & \checkmark & & &\\    
    \hline
    Kernel methods & &\checkmark & &\checkmark & \checkmark & & &\\    
    \hline
    Supervised learning methods & &\checkmark & & & \checkmark & & &\\  
    \hline
    Clustering methods &\checkmark &\checkmark & & \checkmark& & & &\\    
    \hline
    Visualization of high-dimensional data & & & & & & & & \checkmark\\
    \hline    
    Dimension reduction & &\checkmark & & & & & & \\
    \hline
  \end{tabular}
  \caption{\texttt{R} packages for TDA}
  \label{tab:TDA_packages}
\end{table}

The \texttt{TDA} package is the largest \texttt{R} package for TDA. The \texttt{TDA} package offers tools to compute PDs for commonly used types of filtrations such as \emph{Vietoris-Rips}, \emph{Alpha} and \emph{Alpha shape}. It also allows to construct more general sublevel set filtrations and compute the corresponding PDs. Moreover, the \texttt{TDA} package provides implementations to plot PDs and compute \emph{bottleneck and Wasserstein} distances between them. \texttt{TDAstats} offers a variety of tools for conducting statistical inference (such as hypothesis testing) on PDs. Compared to the \texttt{TDA} package, it allows faster computations of PDs for Vietoris-Rips filtrations based on the Ripser C++ library [@Bauer2021Ripser] and more aesthetic visualization of the diagrams using the \texttt{ggplot2} package [@ggplot2]. The \texttt{kernelTDA} package contains implementations of popular kernel-based methods for TDA such as \emph{geodesic Gaussian kernel}, \emph{geodesic Laplacian kernel}, \emph{persistence Fisher kernel} and \emph{persistence sliced Wasserstein kernel}. For computing the Wasserstein distance between a pair of PDs, unlike the \texttt{TDA} package, it uses an iterative procedure to reasonably approximate the exact distance which that leads to a considerable reduction in run-time cost. The \texttt{ripserr} package allows a fast computation of PDs for filtrations on Vietoris-Rips and cubical complexes using the Ripser C++ library. The \texttt{TDApplied} and \texttt{TDAkit} packages offer various tools to integrate topological features (PDs or their vector summaries) into machine and statistical learning settings. The \texttt{tdaunif} is a useful package if one needs to sample points from various manifolds such as a klein-bottle, an ellipse or a torus. \texttt{TDAmapper} offers tools to visualize high-dimensional data by constructing the so-called Mapper graphs that preserve its topological structure.

# Statement of need

Among several vector summaries of PDs available in the TDA literature, only \emph{persistence landscapes}, \emph{persistence silhouettes} [@chazal2014stochastic] and \emph{persistence images} [@adams2017persistence] are implemented in \texttt{R} packages. Moreover, all the code behind the vector implementations is written using standard functions of \texttt{R} which may prove inefficient for large-scale computations. The \texttt{TDAvec} package [@TDAvec] aims to fill in some of these gaps. Its contributions can be summarized in the following three areas: 

1. It expands the list of implemented vector summaries of PDs by providing vectorizations of eight functional summaries found in the TDA literature: \emph{Betti function}\footnote{also called \emph{Betti curve}} [@chazal2021introduction], \emph{persistence landscape function}, \emph{persistence silhouette function}, \emph{persistent entropy summary function} [@atienza2020stability], \emph{Euler characteristic curve} [@richardson2014efficient], \emph{normalized life curve} [@chung2022persistence], \emph{persistence surface} [@adams2017persistence] and \emph{persistence block} [@chan2022computationally].
2. A univariate summary function $f$ of a PD is typically vectorized by evaluating it at each point of a superimposed one dimensional grid and arranging the resulting values into a vector:
\begin{equation}\label{stand_vec}
		(f(t_1),f(t_2),\ldots,f(t_n))\in\mathbb{R}^n,
\end{equation}
where $t_1,t_2,\ldots,t_n$ form an increasing sequence of scale values. For example, the \texttt{landscape()} and \texttt{silhouette()} functions of the \texttt{TDA} package compute vector summaries of persistence landscapes and silhouettes in this manner. The \texttt{TDAvec} package instead employs a different vectorization scheme which involves computing the average values of $f$ between two consecutive scale values $t_i$ and $t_{i+1}$ using integration: 
\begin{equation} 
	\Big(\frac{1}{\Delta t_1}\int_{t_1}^{t_2}f(t)dt,\frac{1}{\Delta t_2}\int_{t_2}^{t_3}f(t)dt,\ldots,\frac{1}{\Delta t_{n-1}}\int_{t_{n-1}}^{t_n}f(t)dt\Big)\in\mathbb{R}^{n-1}, 
\end{equation}
where $\Delta t_i=t_{i+1}-t_i$. Unlike (\ref{stand_vec}), this vectorization method does not miss the behavior of $f$ between neighboring scale points and applies to all univariate summary functions which are easy to integrate, namely persistence silhouette, persistent entropy summary function, Euler characteristic curve, normalized life curve and Betti function. 
3. To achieve higher computational efficiency, all code behind the vector summaries of the \texttt{TDAvec} package is written in C++ using the \texttt{Rcpp} package [@Rcpp]. For example, computing the persistence landscape from a PD with the \texttt{TDAvec} package is more than 200 times faster than with the \texttt{TDA} package.

The \texttt{TDAvec} package and a vignette showing its basic usage with examples are available on the CRAN repository\footnote{https://cran.r-project.org/web/packages/TDAvec/index.html}. 

# Acknowledgements

The author would like to thank Aleksei Luckinsky for his technical assistance to launch the \texttt{TDAvec} package on the CRAN repository.

# References
