import numpy as np
import ripser

from  .tdavec_core import pmin, pmax, DiagToPD, \
    computePersistenceBlock, computePersistenceLandscape, computePersistenceSilhouette, computeNormalizedLife, computeBettiCurve, computeEulerCharacteristic, computePersistentEntropy, computePersistenceImage,\
    computePersistenceBlock_dim0, computePersistenceBlock_dim1, computeFDA, computeAlgebraicFunctions, computeStats, \
    computeComplexPolynomial, computeTemplateFunction, computeTropicalCoordinates

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


def createEllipse(n = 100, a = 1, b = 1, eps = 0.1):
    phi = np.random.uniform(0,2*np.pi,n)
    r = np.random.uniform(1-eps,1+eps,n)
    x = a * r * np.cos(phi)
    y = b * r * np.sin(phi)
    return np.vstack([x,y]).T


class TDAvectorizer:

    def __init__(self, params = {"output": "vpb", "threshold": 2, "inf": "threshold", "maxDim": 1,
                                 "scale": None, "nGrid": 11, "quantiles": False, "tau": 0.3,
                                   "k":1, "sigma": 1, "kFDA":10}):
        self.diags = []
        self.params = params

        self.vectorization_names = [
            "vab", "vpb", "pl", "ps", "nl", "ecc", "pes", "pi", "fda", "algebra", "stats", "poly", 
            "func", "coords"
        ]
        return 

    def getParams(self, varName = None):
        if varName is not None:
            return self.params[varName]
        return self.params
    def setParams(self, params):
        for k in params.keys():
            self.params[k] = params[k]
        return

    def plot(self):
        print("Plotting")
    
    def __getitem__(self, index):
        self.diags = self.diags[index]
    
    def fit(self, data):
        self.diags = []
        for d in data:
            diag_ = ripser.ripser(d, thresh=self.params["threshold"])["dgms"]
            if self.params["inf"] is None:
                diag_[0] = diag_[0][:-1,:]
            elif self.params["inf"] == "threshold":
                diag_[0][-1, 1] = self.params["threshold"]
            else:
                diag_[0][-1, 1] = self.params["inf"]
            self.diags.append(diag_)
        if self.params["scale"] is None:
            limits1 = self.findLimits(1)
            self.params["scale"] = np.linspace(limits1["minB"], limits1["maxB"], self.params["nGrid"])
        return
    
    def findLimits(self, homDim = 0):
        if self.diags is not None:
            births = np.concatenate([d[homDim][:,0] for d in self.diags])
            deaths = np.concatenate([d[homDim][:,1] for d in self.diags])
            return {
                "minB": births.min(), "maxB": births.max(), "minD": deaths.min(), "maxD": deaths.max()
            }

    
    def transform(self, homDim = 1, xSeq = None, ySeq = None, output = None, tau = None, k = None, sigma = None):
        if output is None:
            output = self.params["output"].lower()
        else:
            output = output.lower()
        if xSeq is None:
            xSeq = self.params["scale"]
        if ySeq is None:
            ySeq = xSeq
        if tau is None:
            tau = self.params["tau"]
        if k is None:
            k = self.params["k"]
        if sigma is None:
            sigma = self.params["sigma"]

#    computePersistentEntropy, computePersistenceImage,\
        if type(homDim) == int:
            if output == "vab":
                return np.array([computeBettiCurve(d, homDim = homDim, scaleSeq = xSeq) for d in self.diags])
            elif output == "vpb":
                out = np.array([computePersistenceBlock(d, homDim = homDim, xSeq = xSeq, ySeq = ySeq, tau = tau) for d in self.diags])
                if homDim == 1:
                    out = out.reshape( (out.shape[0], out.shape[1]*out.shape[2]))
                return out
            elif output == "pl":
                return np.array([computePersistenceLandscape(d, homDim = homDim, scaleSeq = xSeq, k = k) for d in self.diags])
            elif output == "ps":
                return np.array([computePersistenceSilhouette(d, homDim = homDim, scaleSeq = xSeq, p=k) for d in self.diags])
            elif output == "nl":
                return np.array([computeNormalizedLife(d, homDim = homDim, scaleSeq = xSeq) for d in self.diags])
            elif output == "ecc":
                return np.array([computeEulerCharacteristic(d, homDim, xSeq) for d in self.diags])
            elif output == "pes":
                return np.array([computePersistentEntropy(d, homDim = homDim, scaleSeq = xSeq) for d in self.diags])
            elif output == "pi":
                return np.array([computePersistenceImage(d, homDim = homDim, xSeq = xSeq, ySeq = ySeq, sigma = sigma) for d in self.diags])
            elif output == "fda":
                maxD = self.findLimits(homDim)["maxD"]
                return np.array([computeFDA(d, maxD, homDim=homDim, K=self.params["kFDA"]) for d in self.diags])
            elif output == "algebra":
                maxD = self.findLimits(homDim)["maxD"]
                return np.array([computeAlgebraicFunctions(d, maxD, homDim=homDim) for d in self.diags])
            elif output == "stats":
                return np.array([computeStats(d, homDim)  for d in self.diags])
            elif output == "poly":
                out = np.array([computeComplexPolynomial(d, homDim)  for d in self.diags])
                out = out[:,0,:]
                return out
            elif output == "func":
                return np.array([computeTemplateFunction(d, homDim)  for d in self.diags])
            elif output == "coords":
                return np.array([computeTropicalCoordinates(d, homDim)  for d in self.diags])

        elif type(homDim) == list:
            out = np.zeros( (len(self.diags), 0) )
            for d in homDim:
                out  = np.hstack([out, 
                                  self.transform(homDim = d, xSeq = xSeq, ySeq = ySeq, output = output, tau = tau, k = k, sigma = sigma)])
            return out
    

def test_package():
    clouds = []
    ratList = np.random.uniform(-0.5, 0.5, 10)
    print("Creating point cloud...")
    for ratio in ratList:
        clouds = clouds + [createEllipse(a=1-ratio, b=1, eps=0.1)]

    vect = TDAvectorizer()
    print("Creating PD...")
    vect.fit(clouds)
    print("Creating vetorization...")
    ps = vect.transform(output="PS", homDim=1)
    return clouds[0], vect.diags[0], ps[0]


    

    

