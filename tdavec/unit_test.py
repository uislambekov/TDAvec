import unittest

from tdavec import TDAvectorizer, tdavec_core, createEllipse
from tdavec.tdavec_core import computeNormalizedLife, computeBettiCurve, computePersistenceBlock, computePersistenceLandscape, computePersistenceSilhouette, \
    computeEulerCharacteristic, computePersistentEntropy, computePersistenceImage, computeComplexPolynomial, computeFDA, computeTropicalCoordinates, computeTemplateFunction
import ripser
import numpy as np

def lists_are_equal(nums1, nums2, atol=1e-6):
    return np.allclose(nums1, nums2, atol=atol)

class Testing_Functions(unittest.TestCase):

    def setUp(self) -> None:
        self.X = np.loadtxt("./tdavec/unitCircle.csv", skiprows=1, delimiter=",")
        self.D = ripser.ripser(self.X, thresh=2)["dgms"]
        self.D[0][-1, 1] = 2
        self.scaleSeq = np.linspace(0, 2, 11)


    def test_PL_0(self):
        python = computePersistenceLandscape(self.D, 0, self.scaleSeq)
        R = [ 0, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2, 0]
        self.assertTrue( lists_are_equal(python, R))
    
    def test_PL_1(self):
        python = computePersistenceLandscape(self.D, 1, self.scaleSeq)
        R = [ 0, 0.0142176303784579, 0.000931983027759931, 0.191114535909928, 
             0.391114535909928, 0.269212150161859, 0.0692121501618592, 0, 0, 0, 0]
        self.assertTrue( lists_are_equal(python, R))

    def test_PS_0(self):
        python = computePersistenceSilhouette(self.D, 0, self.scaleSeq)
        R = [ 0.0507014405880279, 0.0420473568616778, 0.0615450397748265, 0.0861630556847571, 
             0.110781071594688, 0.110781071594688, 0.0861630556847571, 
             0.0615450397748265, 0.0369270238648959, 0.0123090079549653]
        self.assertTrue( lists_are_equal(python, R))

    def test_PS_1(self):
        python = computePersistenceSilhouette(self.D, 1, self.scaleSeq)
        R = [ 0.000119825669467705, 0.00125047122181634, 0.0651147889829073, 0.207585992499712,
              0.257838809390676, 0.120660660329335, 0.00853962588653096, 0, 0, 0]
        self.assertTrue( lists_are_equal(python, R))

    def test_NL_0(self):
        python = computeNormalizedLife(self.D, 0, self.scaleSeq)
        R = [ 0.817130850366702, 0.234100823784605, 0.123090079549653, 0.123090079549653,
              0.123090079549653, 0.123090079549653, 0.123090079549653, 0.123090079549653, 
              0.123090079549653, 0.123090079549653]
        self.assertTrue( lists_are_equal(python, R))

    def test_NL_1(self):
        python = computeNormalizedLife(self.D, 1, self.scaleSeq)
        R = [ 0.0130993953929026, 0.0631855672482697, 0.682417578522476, 
             0.713073264620286, 0.713073264620286, 0.713073264620286, 0.246766669336532, 
             0, 0, 0]
        self.assertTrue( lists_are_equal(python, R))

    def test_VAB_0(self):
        python = computeBettiCurve(self.D, 0, self.scaleSeq)
        R = [ 65.9281367702083, 7.31317883042165, 1, 1, 1, 1, 1, 1, 1, 1]
        self.assertTrue( lists_are_equal(python, R))

    def test_VAB_1(self):
        python = computeBettiCurve(self.D, 1, self.scaleSeq)
        R = [ 0.284203374082668, 1.37424398158854, 1.02801848480822, 1, 1, 1, 0.346060750809296, 0, 0, 0]
        self.assertTrue( lists_are_equal(python, R))

    def test_ECC_0(self):
        python = computeEulerCharacteristic(self.D, 0, self.scaleSeq)
        R = [ 65.9281367702083, 7.31317883042165, 1, 1, 1, 1, 1, 1, 1, 1]
        self.assertTrue( lists_are_equal(python, R))

    def test_ECC_1(self):
        python = computeEulerCharacteristic(self.D, 1, self.scaleSeq)
        R = [ 65.6439333961256, 5.93893484883312, -0.0280184848082246, 0, 0, 0, 0.653939249190704, 1, 1, 1]
        self.assertTrue( lists_are_equal(python, R))

    def test_PES_0(self):
        python = computePersistentEntropy(self.D, 0, self.scaleSeq)
        R = [ 4.82683006650221, 1.01731577228001, 0.372004512741568, 0.372004512741568, 0.372004512741568, 
             0.372004512741568, 0.372004512741568, 0.372004512741568, 0.372004512741568, 0.372004512741568]
        self.assertTrue( lists_are_equal(python, R))

    def test_PES_1(self):
        python = computePersistentEntropy(self.D, 1, self.scaleSeq)
        R = [ 0.0567767350968682, 0.264582251401072, 0.338458609058045, 0.347892602095769, 0.347892602095769, 
             0.347892602095769, 0.120391975082262, 0, 0, 0]
        self.assertTrue( lists_are_equal(python, R))

    def test_FDA_0(self):
        maxD = max(self.D[0][:,1])
        python = computeFDA(self.D, maxD, 0, 10)
        R = [ 8.12413156006299, 2.78355597654528, 9.43745477301271, 4.90077947666511, 
             7.75038019279065, 5.9773445200531, 5.57180752169283, 6.05176570843183, 3.51266625324211, 
             5.47299594705235, 1.97258501571194, 4.66807567491559, 1.03655377773163, 3.93750799080157, 
             0.554577049981595, 3.38441743756259, 0.309407953491623, 2.97285589648734, 0.148884628338501, 
             2.63191282881754, 0.0204442278795615]
        self.assertTrue( lists_are_equal(python, R))

    def test_FDA_1(self):
        maxD = max(self.D[1][:,1])
        python = computeFDA(self.D, maxD, 1, 10)
        R = [ 0.950593892521343, 0.0293671377520018, -0.133034606398945, -0.0675654244183837, -0.178134342570056,
              -0.158710242349944, -0.105821960064297, -0.163135297451828, 0.00444685439650682, -0.0813024955664044,
                0.0748359974671063, -0.00265992760211067, 0.0341593046222416, -0.0289683485110909, -0.00577124724551058,
                  -0.0232309859037072, 0.0331392480951682, 0.0289202994258002, -0.000687313921102223, -0.0132191587942031, 
                  -0.0499634291792129]
        self.assertTrue( lists_are_equal(python, R))

    def test_VPB_0(self):
        ySeqH0 = np.quantile(self.D[0][:,1] - self.D[0][:,0] , np.arange(0, 1.1, 0.2))
        vpb0 = computePersistenceBlock(self.D, homDim=0, xSeq = [], ySeq = ySeqH0)
        R = [ 0.809001048929071, 3.64481540058609, 5.48661206724003, 6.38281631026758, 1.01395250734014]
        self.assertTrue( lists_are_equal(vpb0, R))

    def test_VPB_1(self):
        xSeqH1 = np.quantile(self.D[1][:,0], np.arange(0, 1.1, 0.2))
        ySeqH1 = np.quantile(self.D[1][:,1]- self.D[1][:,0], np.arange(0, 1.1, 0.2))
        vpb1 = computePersistenceBlock(self.D, homDim = 1, xSeq=xSeqH1, ySeq=ySeqH1)
        vpb1 = np.transpose(vpb1).reshape( (25,))
        R = [ 0, 0.000499559782281428, 0.00903377330326861, 0.0305098817964819, 
             0.00579736873795015, 0.00951944957273334, 0.0261053156675721, 0.0203673036116531, 
             0.00428701302896287, 0.0107446149362357, 0.0425768331028477, 0.105346832126612,
               0.145880724120365, 0, 0, 0.0328035481022359, 0.0276785592055745, 0.108908859370022, 
               0.131882225448411, 0.0168304894683267, 0.293939456124958, 0.316298612637143, 0.334451098933436,
                 0.356748693217673, 0.36362269411629] 
        self.assertTrue( lists_are_equal(vpb1, R))

    def test_PI_0(self):
        resB, resP = 5, 5
        minPH0, maxPH0 = np.min(self.D[0][:,1]), np.max(self.D[0][:,1])
        ySeqH0 = np.linspace(minPH0, maxPH0, resP+1)
        xSeqH0 = np.zeros( resB+1)
        sigma = 0.5*(maxPH0-minPH0)/resP
        pi0 = computePersistenceImage(self.D, homDim = 0, xSeq = xSeqH0, ySeq = ySeqH0, sigma = sigma)
        pi0_R = [ 4.56670703042467, 0.965627354019496, 0.0113500659339766, 0.0227233996700507, 0.477249868112126]
        self.assertTrue( lists_are_equal(pi0, pi0_R))

    def test_PI_1(self):
        resB, resP = 5, 5
        minBH1, maxBH1 = np.min(self.D[1][:,0]), np.max(self.D[1][:,0])
        xSeqH1 = np.linspace(minBH1, maxBH1, resB+1)
        minPH1, maxPH1 = np.min(self.D[1][:,1] - self.D[1][:,0]), np.max(self.D[1][:,1] - self.D[1][:,0])
        ySeqH1 = np.linspace(minPH1, maxPH1, resP+1)
        sigma = 0.5*(maxPH1-minPH1)/resP
        pi1 = computePersistenceImage(self.D, homDim = 1, xSeq = xSeqH1, ySeq = ySeqH1, sigma = sigma)
        pi1_R = [ 0.0437160130945874, 0.0495556729840402, 0.0451237783565141, 0.0329995660586893, 0.0186444239574634, 
                 0.00752765376757397, 0.00776084705021553, 0.00630598644271661, 0.00425097166237492, 0.00227854726438653, 
                 6.03800613955066e-05, 5.84101446147877e-05, 4.42739912785808e-05, 3.09472597878289e-05, 
                 2.02106358360438e-05, 0.000183467522328633, 0.00089210088336952, 0.00271421546306594, 
                 0.00517182299303047, 0.00617537792830458, 0.00385383980210524, 0.0187402270298275, 0.0570177372523999, 
                 0.108645123651971, 0.12972698923131] 
        self.assertTrue( lists_are_equal(pi1, pi1_R))

    def test_complex_poly_0(self):
        poly0 = computeComplexPolynomial(self.D, homDim = 0)
        poly0_R = np.array([0, -16.248263120126])
        self.assertTrue( lists_are_equal(poly0, poly0_R))

    def test_complex_poly_1(self):
        poly1 = computeComplexPolynomial(self.D, homDim = 1)
        poly1_R = np.array([-3.81072113328258, -5.01722645154033])
        self.assertTrue( lists_are_equal(poly1, poly1_R))

    def test_tropical_coords_0(self):
        dat = computeTropicalCoordinates(self.D, homDim = 0)
        dat_R = np.array([2, 2.36289119595766, 2.68205181705815, 2.98638484543615, 16.248263120126, 0, 183.751736879874])
        self.assertTrue( lists_are_equal(dat[:-1], dat_R[:-1]))
    
    def test_tropical_coords_1(self):
        dat = computeTropicalCoordinates(self.D, homDim = 1)
        dat_R = np.array([ 0.860326686071788, 0.954869453646198, 1.02284830582612, 1.0654713502642, 1.20650531825775, 0.755064096276033, 14.5381885375704])
        self.assertTrue( lists_are_equal(dat[:-1], dat_R[:-1]))

    def test_temp_func_0(self):
        dat = computeTemplateFunction(self.D, homDim = 0)
        dat_R = np.array( [47.7068836896669, 27.5322316257831, 9.10973719256688, 0.620518170581516, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0999999999999992, 0.900000000000002])
        self.assertTrue( lists_are_equal(dat[:-1], dat_R[:-1]))

    def test_temp_func_1(self):
        dat = computeTemplateFunction(self.D, homDim = 1)
        dat_R = np.array([
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.468406748165331, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.75451666854592, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13899057979361, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.365103088717902, 0, 0, 0, 0, 0, 0, 0.496733139282122, 0.503266860717878, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0355717802395706, 0, 0, 0, 0, 0, 0, 0.0888546409007165, 0.0888546409007165, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
        ])
        self.assertTrue( lists_are_equal(dat[:-1], dat_R[:-1]))

class Testing_Class(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.clouds = []
        self.ratList = np.random.uniform(-0.5, 0.5, 10)
        for ratio in self.ratList:
            self.clouds = self.clouds + [createEllipse(a=1-ratio, b=1, eps=0.1)]
        self.vect = TDAvectorizer()
        self.vect.fit(self.clouds)

    def test_N(self):
        self.assertEqual( len(self.vect.diags), 10)

    def test_dim_PS(self):
        xSeq = np.linspace(-1, 1, 21)
        ps = self.vect.transform(output = "ps", homDim = 0, xSeq = xSeq)
        self.assertEqual( ps.shape, (10, 20))




if __name__ == '__main__':
    unittest.main()