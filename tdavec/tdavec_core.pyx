#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
cimport numpy as np
from libc.math cimport log2
from cython cimport boundscheck, wraparound, nonecheck
from libc.math cimport sqrt, cos, sin
from libc.math cimport fabs, fmax

def DiagToPD(D):
    """
    Generates a list of persistence diagrams (PD) from a given list of persistence diagrams (D).

    Parameters:
    - D (list): A list of persistence diagrams, where each diagram is represented as a numpy array.

    Returns:
    - PD (list): A list of persistence diagrams (PD), where each PD is represented as a numpy array.
      Each PD contains two columns: the first column represents the birth values of the persistence pairs,
      and the second column represents the persistence, i.e. the death values minus the birth values.
    """
    PD = [ np.transpose(np.array([D[dim][:,0], D[dim][:,1] - D[dim][:,0]])) for dim in range(len(D))]
    return PD

def dict_to_array(d, to_list = True):
    if to_list:
        return np.array(list(d.values()))
    else:
        return d

def computePersistenceBlock_dim0(np.ndarray[np.float64_t, ndim=1] x, 
                    np.ndarray[np.float64_t, ndim=1] y, 
                    np.ndarray[np.float64_t, ndim=1] ySeq, 
                    np.ndarray[np.float64_t, ndim=1] lam):
    """
    Compute the VPB values for dimension 0.

    Parameters:
        x (numpy.ndarray): The x values.
        y (numpy.ndarray): The y values.
        ySeq (numpy.ndarray): The sequence of y values.
        lam (numpy.ndarray): The lambda values.

    Returns:
        numpy.ndarray: The computed VPB values.
    """
    cdef Py_ssize_t n = ySeq.shape[0] - 1
    cdef int nPoints = y.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] vpb = np.zeros(n, dtype=np.float64)
    cdef int i, j
    cdef double c, d, y_cd, lam_cd, yMin, yMax
    for i in range(n):
        c = ySeq[i]
        d = ySeq[i+1]
        for j in range(nPoints):
            if c - lam[j] < y[j] < d + lam[j]:
                y_cd = y[j]
                lam_cd = lam[j]
                yMin = max(c, y_cd - lam_cd)
                yMax = min(d, y_cd + lam_cd)
                vpb[i] += 0.5 * (yMax**2 - yMin**2) / (ySeq[i+1] - ySeq[i])
    return vpb

def pmax(num, vec):
    """
    Compute the element-wise maximum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise maximum.
    """
    return np.array([max(num, vec[i_]) for i_ in range(vec.size)])

def pmin(num, vec):
    """
    Compute the element-wise minimum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise minimum.
    """
    return np.array([min(num, vec[i_]) for i_ in range(vec.size)])

def computePersistenceBlock_dim1(np.ndarray[np.float64_t, ndim=1] x, 
                    np.ndarray[np.float64_t, ndim=1] y, 
                    np.ndarray[np.float64_t, ndim=1] xSeq, 
                    np.ndarray[np.float64_t, ndim=1] ySeq, 
                    np.ndarray[np.float64_t, ndim=1] lam):
    """
    Compute the Vector Persistence Block (VPB) vectorization for a given set of points in dimension 1.

    Parameters:
        x (numpy.ndarray): The x-coordinates of the points.
        y (numpy.ndarray): The y-coordinates of the points.
        xSeq (numpy.ndarray): The x-coordinates of the grid points.
        ySeq (numpy.ndarray): The y-coordinates of the grid points.
        lam (numpy.ndarray): The lambda values.

    Returns:
        numpy.ndarray: The VPB matrix.
    """
    cdef int n = xSeq.shape[0] - 1
    cdef int m = ySeq.shape[0] - 1
    cdef np.ndarray[np.float64_t, ndim=1] dx = np.diff(xSeq)
    cdef np.ndarray[np.float64_t, ndim=1] dy = np.diff(ySeq)
    cdef np.ndarray[np.float64_t, ndim=2] vpb = np.zeros((n, m), dtype=np.float64)
    cdef int i, j, k, num_inds
    cdef double a, b, c, d, add
    cdef np.ndarray[np.int64_t, ndim=1] inds
    cdef np.ndarray[np.float64_t, ndim=1] xInd, yInd, lamInd, xMin, xMax, yMin, yMax

    for i in range(n):
        a, b = xSeq[i], xSeq[i+1]
        for j in range(m):
            c, d = ySeq[j], ySeq[j+1]
            # Using bitwise operations for logical conditions
            xCond = (x + lam >= a) & (x - lam <= b)
            yCond = (y + lam >= c) & (y - lam <= d)
            inds = np.where(xCond & yCond)[0]
            num_inds = inds.shape[0]

            if num_inds > 0:
                xInd = x.take(inds)
                yInd = y.take(inds)
                lamInd = lam.take(inds)
                xMin = np.maximum(a, xInd - lamInd)
                xMax = np.minimum(b, xInd + lamInd)
                yMin = np.maximum(c, yInd - lamInd)
                yMax = np.minimum(d, yInd + lamInd)
                add = 0.5 * np.sum((xMax - xMin) * (yMax - yMin) * (xMax + xMin + yMax + yMin)) / dx[i] / dy[j]
                vpb[i, j] += add
    return np.asarray(vpb)

def computePersistenceBlock(D, homDim, xSeq, ySeq, tau=0.3):
    """
    Compute the VPB vectorization using the given parameters.

    Parameters:
        D (list): Persistence Diagram (list of birth-death arrays for each dimension).
        homDim (int): The dimension along which the homogeneity is computed.
        xSeq (numpy.ndarray): The x-coordinates of the grid points.
        ySeq (numpy.ndarray): The y-coordinates of the grid points.
        tau (float, optional): The tau value. Defaults to 0.3.

    Returns:
        numpy.ndarray: The VPB matrix.
    """
    x = D[homDim][:,0]
    y = D[homDim][:,1] - x
    lam = tau * y
    if homDim == 0:
        return computePersistenceBlock_dim0(x, y, ySeq, lam)
    else:
        return computePersistenceBlock_dim1(x, y, xSeq, ySeq, lam)

def computePersistenceLandscape(D, homDim, scaleSeq, k=1):
    """
    Compute the persistence landscape (PL) for a given homological dimension, scale sequence, and order of landscape.

    Parameters:
        D (numpy.ndarray): Persistence Diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the PL is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.
        k (int, optional): The order of the PL. Defaults to 1.

    Returns:
        numpy.ndarray: The persistence landscape vector.
    """
    if D[homDim].shape[0] == 0:
        return np.zeros( len(scaleSeq))
    birth, death = D[homDim][:,0], D[homDim][:,1]
    Lambda = [
        np.sort(pmax(0, np.apply_along_axis(min, 0, np.array([s-birth, death-s]))))[-k]
        for s in scaleSeq]
    return np.array(Lambda)

def computePersistenceSilhouette(D, homDim, scaleSeq, p=1):
    """
    Compute the Persistence Silhouette vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the PS is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.
        p (int, optional): The power to raise the difference between y and x. Defaults to 1.

    Returns:
        numpy.ndarray: The persistence spectrum vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    pp = (y-x)**p
    w = pp/np.sum(pp)

    phi = []
    for k in range(len(scaleSeq)-1):
        alpha1 = pmax(scaleSeq[k], x)
        alpha2 = pmax(scaleSeq[k], (x+y)/2)
        beta1 = pmin(scaleSeq[k+1], (x+y)/2)
        beta2 = pmin(scaleSeq[k+1], y)
        b1 = pmax(0,beta1-alpha1)*((beta1+alpha1)/2-x)
        b2 = pmax(0,beta2-alpha2)*(y-(beta2+alpha2)/2)
        phi.append( np.sum(w*(b1+b2))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(phi)

def computeNormalizedLife(D, homDim, scaleSeq):
    """
    Compute theNormalized Life Curve vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the NL is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The nonlinear landscape vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    lL = (y-x)/sum(y-x)
    nl = []
    for k in range(len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        nl.append( np.sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(nl)

def computeBettiCurve(D, homDim, scaleSeq):
    """
    Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the VAB is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The VAB vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    vab = []
    for k in range( len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(vab)

def computeEulerCharacteristic(D, maxhomDim, scaleSeq):
    """
    Compute the Euler Characteristic Curve (ECC) vectorization for a given homological dimension, maximum homological dimension, and scale sequence.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        maxhomDim (int): The maximum homological dimension.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The ECC vector.
    """
    ecc = np.zeros( len(scaleSeq)-1)
    for d in range(maxhomDim+1):
        ecc = ecc + (-1)**d * computeBettiCurve(D, d, scaleSeq)
    return ecc

def computePersistentEntropy(D, homDim, scaleSeq):
    """
    Compute the Persistence Entropy Summary (PES) vectorization for a given homological dimension, scale sequence, and persistence diagram.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        list: The PES values.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    lL = (y-x)/np.sum(y-x)
    entr = -lL*np.log10(lL)/np.log10(2)
    pes = []
    for k in range( len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        pes.append( np.sum(entr*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return pes

from scipy.stats import norm
def pnorm(x, mean, sd):
    """
    Calculate the cumulative distribution function of a normal distribution.

    Parameters:
        x (float): The value at which to calculate the cumulative distribution function.
        mean (float): The mean of the normal distribution.
        sd (float): The standard deviation of the normal distribution.

    Returns:
        float: The cumulative distribution function value at x.
    """
    return norm.cdf(x, mean, sd)

def outer(x, y):
    """
    Generate the outer product of two arrays.

    Parameters:
        x (array-like): The first input array.
        y (array-like): The second input array.

    Returns:
        numpy.ndarray: The outer product of the input arrays.
    """
    return np.array([x_*y_ for y_ in y for x_ in x])

def PSurfaceH0(point, y_lower, y_upper, sigma, maxP):
    """
    Calculate the surface probability density function for a specific homDim=0 point on the y-axis .

    Parameters:
        point (tuple): A tuple containing the x and y coordinates of the point.
        y_lower (float): The lower bound of the y-axis interval.
        y_upper (float): The upper bound of the y-axis interval.
        sigma (float): The standard deviation of the normal distribution.
        maxP (float): The maximum value of the y-axis.

    Returns:
        float: The surface probability density function value at the given point.

    """
    y = point[1]
    out2 = pnorm(y_upper, y, sigma) - pnorm(y_lower, y, sigma)
    wgt = y/maxP if y<maxP else 1
    return wgt*out2

def PSurfaceHk(point, y_lower, y_upper, x_lower, x_upper, sigma, maxP):
    """
    Calculate the surface probability density function for a specific homDim>0point in a two-dimensional space.

    Parameters:
        point (tuple): A tuple containing the x and y coordinates of the point.
        y_lower (float): The lower bound of the y-axis interval.
        y_upper (float): The upper bound of the y-axis interval.
        x_lower (float): The lower bound of the x-axis interval.
        x_upper (float): The upper bound of the x-axis interval.
        sigma (float): The standard deviation of the normal distribution.
        maxP (float): The maximum value of the y-axis.

    Returns:
        float: The surface probability density function value at the given point.

    """
    x, y = point[0], point[1]
    out1 = pnorm(x_upper,x,sigma) - pnorm(x_lower,x,sigma)
    out2 = pnorm(y_upper,y,sigma) - pnorm(y_lower,y,sigma)
    wgt = y/maxP if y<maxP else 1
    return wgt*outer(out1, out2)

def computePersistenceImage(D, homDim, xSeq, ySeq, sigma):
    """
    Compute the surface Persistence Image (PI) vectorization for a given Persistence diagram

    Args:
        D (list): Persistence Diagram (list of birth-death arrays for each dimension).
        homDim (int): The dimension to compute the surface probability density function for.
        xSeq (list): The x-axis sequence.
        ySeq (list): The y-axis sequence.
        sigma (float): The standard deviation of the normal distribution.

    Returns:
        numpy.ndarray: The surface probability density function values for each data point.

    """
    D_ = D[homDim].copy()
    D_[:,1] = D_[:,1] - D_[:,0]
    n_rows = D_.shape[0]

    resB = len(xSeq) - 1
    resP = len(ySeq)-1
    minP, maxP = ySeq[0], ySeq[-1]
    dy = (maxP-minP)/resP
    y_lower = np.arange(minP, maxP, dy)
    y_upper = y_lower + dy

    nSize = resP if homDim == 0 else resP*resB
    Psurf_mat = np.zeros( (nSize, n_rows))
    if homDim==0:
        for i in range(n_rows):
            Psurf_mat[:, i] = PSurfaceH0(D_[i, :], y_lower, y_upper, sigma, maxP)
    else:
        minB, maxB = xSeq[0], xSeq[-1]
        dx = (maxB-minB)/resB
        x_lower = np.arange(minB, maxB, dx)
        x_upper = x_lower + dx
        for i in range(n_rows):
            Psurf_mat[:, i] = PSurfaceHk(D_[i, :], y_lower, y_upper, x_lower, x_upper, sigma, maxP)
    out = np.sum(Psurf_mat, axis = 1)
    return out

def computeFDA(PD, maxD, homDim = 0, K = 10):
    X = np.zeros( (2*K+1))
    pd = PD[homDim]
    b = pd[:,0]/maxD; d = pd[:,1]/maxD
    X[0] = np.sum(d - b)
    for m in range(1, K+1):
        c = 2*m*np.pi
        alpha_sin = np.sin(c*d)-np.sin(c*b)
        alpha_cos = np.cos(c*d)-np.cos(c*b)
        X[2*m-1] = -np.sqrt(2)/c * np.sum(alpha_cos)
        X[2*m] = np.sqrt(2)/c * np.sum(alpha_sin)
    return X

def computeAlgebraicFunctions(PD, maxD, homDim = 0):
    pd = PD[homDim]
    pers = pd[:,1]-pd[:,0]
    return(np.array([
        np.sum(pd[:,0]*pers),
        np.sum( (maxD - pd[:,1])*pers),
        np.sum( pd[:,0]**2 * pers**4),
        np.sum( (maxD - pd[:,1])**2*pers**4)
    ]))


#========= computeStats
# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False

# Expose this to Python
def computeStats(list diagrams, int hom_dim, to_list = True):
    if hom_dim < 0 or hom_dim >= len(diagrams):
        raise ValueError("Invalid homological dimension index.")

    cdef np.ndarray[np.float64_t, ndim=2] D = diagrams[hom_dim]
    cdef Py_ssize_t n_rows = D.shape[0]

    if n_rows == 0:
        return _empty_result(to_list)

    cdef np.float64_t[:, :] D_view = D
    cdef list x_list = []
    cdef list y_list = []
    cdef int i
    cdef int j
    cdef double val


    for i in range(n_rows):
        if np.isfinite(D_view[i, 1]):
            x_list.append(D_view[i, 0])
            y_list.append(D_view[i, 1])

    cdef int total_bars = len(x_list)
    if total_bars == 0:
        return _empty_result(to_list)

    cdef np.ndarray[np.float64_t, ndim=1] x = np.array(x_list, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] y = np.array(y_list, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] mid = (x + y) / 2
    cdef np.ndarray[np.float64_t, ndim=1] l = y - x

    cdef stats_births = _calc_stats(x)
    cdef stats_deaths = _calc_stats(y)
    cdef stats_midpoints = _calc_stats(mid)
    cdef stats_lifespans = _calc_stats(l)

    cdef double L = np.sum(l)
    cdef double entropy = 0.0
    if L > 0:
        for j in range(l.shape[0]):
            val = l[j] / L
            entropy -= val * log2(val)

    cdef dict out = {}
    cdef list keys = ['births', 'deaths', 'midpoints', 'lifespans']
    cdef list names = ['mean', 'stddev', 'median', 'iqr', 'range', 'p10', 'p25', 'p75', 'p90']
    cdef int s_idx
    for k, stats in zip(keys, [stats_births, stats_deaths, stats_midpoints, stats_lifespans]):
        for s_idx, name in enumerate(names):
            out[f"{name}_{k}"] = stats[s_idx]

    out["total_bars"] = total_bars
    out["entropy"] = entropy
    return dict_to_array(out, to_list)


cdef np.ndarray[np.float64_t, ndim=1] _calc_stats(np.ndarray[np.float64_t, ndim=1] values):
    cdef np.ndarray[np.float64_t, ndim=1] result = np.zeros(9, dtype=np.float64)
    if values.shape[0] < 2:
        return result

    cdef np.ndarray[np.float64_t, ndim=1] P = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    cdef np.ndarray[np.float64_t, ndim=1] qt = np.quantile(values, P)

    result[0] = np.mean(values)
    result[1] = np.std(values, ddof=1)
    result[2] = qt[3]  # Median
    result[3] = qt[4] - qt[2]  # IQR
    result[4] = qt[6] - qt[0]  # Range
    result[5] = qt[1]  # 10th percentile
    result[6] = qt[2]  # 25th percentile
    result[7] = qt[4]  # 75th percentile
    result[8] = qt[5]  # 90th percentile

    return result


cdef _empty_result(to_list = True):
    cdef dict d = {}
    cdef list names = ['mean', 'stddev', 'median', 'iqr', 'range', 'p10', 'p25', 'p75', 'p90']
    cdef list keys = ['births', 'deaths', 'midpoints', 'lifespans']
    for k in keys:
        for n in names:
            d[f"{n}_{k}"] = 0.0
    d["total_bars"] = 0
    d["entropy"] = 0.0
    return dict_to_array(d, to_list)



# Allow complex128 arrays
ctypedef np.complex128_t complex_t

# Helper function S(x, y)
def S(np.ndarray[np.float64_t, ndim=1] x,
      np.ndarray[np.float64_t, ndim=1] y):
    cdef int n = x.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] alpha = np.sqrt(x**2 + y**2)
    cdef np.ndarray[np.float64_t, ndim=1] factor = (y - x) / (alpha * np.sqrt(2.0))
    factor[~np.isfinite(alpha)] = 0.0
    return factor * x + 1j * (factor * y)

def T(np.ndarray[np.float64_t, ndim=1] x,
      np.ndarray[np.float64_t, ndim=1] y):
    cdef int n = x.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] alpha = np.sqrt(x**2 + y**2)
    cdef np.ndarray[np.float64_t, ndim=1] cos_alpha = np.cos(alpha)
    cdef np.ndarray[np.float64_t, ndim=1] sin_alpha = np.sin(alpha)
    cdef np.ndarray[np.float64_t, ndim=1] factor = (y - x) / 2.0
    return factor * (cos_alpha - sin_alpha) + 1j * factor * (cos_alpha + sin_alpha)
# Main function
def computeComplexPolynomial(list diag,
                             int homDim,
                             int m=1,
                             str polyType="R"):

    if homDim >= len(diag):
            raise ValueError("homDim exceeds number of diagram dimensions")

    cdef np.ndarray[double, ndim=2] D = diag[homDim]
    if D.shape[1] != 2:
        raise ValueError("Each diagram must be an (n, 2) array")

    x = D[:, 0]
    y = D[:, 1]

    # Remove non-finite death times
    finite_mask = np.isfinite(y)
    x = x[finite_mask]
    y = y[finite_mask]

    if x.shape[0] == 0:
        return np.zeros((m, 2), dtype=np.float64)

    if x.shape[0] < m:
        raise ValueError("m must be less than or equal to the number of points in the diagram!")

    # Compute complex roots
    if polyType == "R":
        roots = x + 1j * y
    elif polyType == "S":
        roots = S(x, y)
    elif polyType == "T":
        roots = T(x, y)
    else:
        raise ValueError("Choose between polyType = 'R', 'S', or 'T'.")

    # Polynomial multiplication
    poly = np.array([1.0 + 0j])  # Start with constant 1
    for root in roots:
        poly = np.convolve(poly, np.array([1.0, -root]))

    # Get m coefficients, skipping the constant term
    real_part = np.real(poly[1:m+1])
    imag_part = np.imag(poly[1:m+1])

    return np.column_stack((real_part, imag_part))

# distutils: language = c++
# cython: boundscheck=False, wraparound=False, cdivision=True

def tent_function_1D(np.ndarray[np.float64_t, ndim=1] y, double b, double delta) -> double:
    cdef Py_ssize_t i, n = y.shape[0]
    cdef double result = 0.0
    cdef double diff
    for i in range(n):
        diff = fabs(y[i] - b)
        result += max(0.0, 1.0 - (diff / delta))
    return result

def tent_function_2D(np.ndarray[np.float64_t, ndim=1] x,
                     np.ndarray[np.float64_t, ndim=1] y,
                     double a, double b, double delta) -> double:
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double result = 0.0
    cdef double dx, dy, max_dist
    for i in range(n):
        dx = fabs(x[i] - a)
        dy = fabs(y[i] - b)
        max_dist = max(dx, dy)
        result += max(0.0, 1.0 - (max_dist / delta))
    return result

def computeTemplateFunction(list diagram, int homDim, double delta = 0.1, int d=20, double epsilon=0.01):
    if delta < 0 or epsilon < 0:
        raise ValueError("The arguments 'delta' and 'epsilon' must be positive!")

    cdef np.ndarray[np.float64_t, ndim=2] D = diagram[homDim]
    cdef np.ndarray[np.float64_t, ndim=1] x = D[:, 0]
    cdef np.ndarray[np.float64_t, ndim=1] y = D[:, 1]
    cdef np.ndarray[np.float64_t, ndim=1] tf = np.zeros(d, dtype=np.float64)

    # Remove entries with non-finite death times
    mask = np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if np.any(x < 0):
        raise ValueError("The birth values must all be positive!")

    cdef np.ndarray[np.float64_t, ndim=1] l = y - x
    cdef Py_ssize_t i, j, idx
    cdef double a, b

    if x.size == 0:
        return np.zeros((d+1)*d, dtype=np.float64)

    sumX = np.sum(np.abs(np.diff(x)))
    if homDim == 0 and sumX == 0:
        center_l = np.linspace(delta, d * delta, d) + epsilon
        for j in range(d):
            tf[j] = tent_function_1D(l, center_l[j], delta)
        return tf

    cdef np.ndarray[np.float64_t, ndim=1] tf2 = np.zeros((d+1)*d, dtype=np.float64)
    center_x = np.linspace(0, d * delta, d+1)
    center_l = np.linspace(delta, d * delta, d) + epsilon

    idx = 0
    for i in range(d+1):
        for j in range(d):
            a = center_x[i]
            b = center_l[j]
            tf2[idx] = tent_function_2D(x, l, a, b, delta)
            idx += 1
    return tf2


@boundscheck(False)
@wraparound(False)
def computeTropicalCoordinates(list D, int homDim, int r=1, to_list = True):
    """
    Parameters:
    - D: list of numpy arrays, one per homological dimension
    - homDim: the homological dimension to use
    - r: scalar multiplier for lifespans
    Returns:
    - dict with keys F1 through F7
    """
    if r <= 0:
        raise ValueError("r must be a positive integer!")

    cdef out  = dict()
    cdef out0 = {f"F{i+1}": 0.0 for i in range(7)}
    out0 = dict_to_array(out0, to_list) 

    if homDim >= len(D):
        return out0

    cdef np.ndarray[np.float64_t, ndim=2] data = D[homDim]

    if data.shape[0] == 0:
        return out0

    # Extract x and y
    cdef np.ndarray[np.float64_t, ndim=1] x = data[:, 0]
    cdef np.ndarray[np.float64_t, ndim=1] y = data[:, 1]

    finite_idx = np.isfinite(y)
    x = x[finite_idx]
    y = y[finite_idx]


    if x.shape[0] == 0:
        out =  out0

    cdef np.ndarray[np.float64_t, ndim=1] lambda_ = y - x
    cdef np.ndarray[np.float64_t, ndim=1] l = np.sort(lambda_)[::-1]

    # F1 to F4
    cdef double F1 = l[0]
    cdef double F2, F3, F4
    cdef Py_ssize_t n = l.shape[0]

    if n > 3:
        F2 = l[0] + l[1]
        F3 = F2 + l[2]
        F4 = F3 + l[3]
    elif n == 3:
        F2 = l[0] + l[1]
        F3 = F2 + l[2]
        F4 = F3
    elif n == 2:
        F2 = l[0] + l[1]
        F3 = F2
        F4 = F2
    else:
        F2 = F1
        F3 = F1
        F4 = F1

    F5 = np.sum(l)
    d = np.minimum(r * lambda_, x)
    F6 = np.sum(d)
    F7 = np.sum(np.maximum(d + lambda_, 0) - (d + lambda_))

    out =  {
        "F1": F1,
        "F2": F2,
        "F3": F3,
        "F4": F4,
        "F5": F5,
        "F6": F6,
        "F7": F7
    }

    out = dict_to_array(out, to_list)

    return out