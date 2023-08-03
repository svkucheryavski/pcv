import numpy as np

import unittest
import itertools
import math
import os

from src.pcv.methods import pcvpca, get_pcamodel, get_xpvorth

def pca_predict(X: np.ndarray, Xpv: np.ndarray, ncomp: int, center: bool = True, scale: bool = False):
    """
    Create a global PCA model using dataset "Xc", applies this model to dataset "Xpv"
    and return matrices with scoures (T), score distance (H) and orthogonal distance (Q).
    """

    (nrows, ncols) = X.shape

    # compute global mean and standard deviation and autoscale the whole data
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    X = (X - mX) / sX

    # create a global model
    P, s = get_pcamodel(X, ncomp)

    # apply the model
    Xpv = (Xpv - mX) / sX
    Tpv = np.dot(Xpv, P)
    Upv = Tpv / (s / math.sqrt(nrows - 1))

    Hpv = np.zeros((nrows, ncomp))
    Qpv = np.zeros((nrows, ncomp))

    for a in range(1, ncomp + 1):
        Pa = P[..., :a]
        Tpva = Tpv[..., :a]
        Upva = Upv[..., :a]
        Epva = Xpv - np.dot(Tpva, np.transpose(Pa))
        Qpv[..., a - 1] = (Epva * Epva).sum(axis = 1)
        Hpv[..., a - 1] = (Upva * Upva).sum(axis = 1)

    return {'T': Tpv, 'H': Hpv, 'Q': Qpv}



def test_pcacase(X, ncomp, cv, center, scale, T, H, Q):
    """
    Create PCA global model, apply it to generated Xpv set and test equality of
    main outcomes (scores, score distances and orthogonal distances) entered manually
    from R tests.
    """

    Xpv = pcvpca(X, ncomp, center = center, scale = scale, cv = cv)
    r = pca_predict(X, Xpv, ncomp, center = center, scale = scale)
    np.testing.assert_array_almost_equal(r['T'], np.array(T).reshape(X.shape[0], ncomp))
    np.testing.assert_array_almost_equal(r['H'], np.array(H).reshape(X.shape[0], ncomp))
    np.testing.assert_array_almost_equal(r['Q'], np.array(Q).reshape(X.shape[0], ncomp))



def test_pcacase_ref(X, ncomp, cv, center, scale):
    """
    Create PCA global model, apply it to generated Xpv set and test equality of
    main outcomes (scores, score distances and orthogonal distances) taken from a
    reference file generated from R tests
    """

    file_suffix = '-' + str(ncomp) + '-' + str(scale).upper() + '-' + \
        (cv['type'] + str(cv['nseg']) if 'nseg' in cv else cv['type']) + '.csv'

    Qg = np.genfromtxt('../.tests/pcvpca/Qpvg' + file_suffix, delimiter=',')
    Hg = np.genfromtxt('../.tests/pcvpca/Hpvg' + file_suffix, delimiter=',')
    Ql = np.genfromtxt('../.tests/pcvpca/Qpvl' + file_suffix, delimiter=',')
    Hl = np.genfromtxt('../.tests/pcvpca/Hpvl' + file_suffix, delimiter=',')

    Xpvg = pcvpca(X, ncomp, center = center, scale = scale, cv = cv, cvscope = 'global')
    Xpvl = pcvpca(X, ncomp, center = center, scale = scale, cv = cv, cvscope = 'local')

    rg = pca_predict(X, Xpvg, ncomp, center = center, scale = scale)
    rl = pca_predict(X, Xpvl, ncomp, center = center, scale = scale)

    np.testing.assert_array_almost_equal(rg['Q'].reshape(X.shape[0], ncomp), Qg.reshape(X.shape[0], ncomp))
    np.testing.assert_array_almost_equal(rg['H'].reshape(X.shape[0], ncomp), Hg.reshape(X.shape[0], ncomp))
    np.testing.assert_array_almost_equal(rl['Q'].reshape(X.shape[0], ncomp), Ql.reshape(X.shape[0], ncomp))
    np.testing.assert_array_almost_equal(rl['H'].reshape(X.shape[0], ncomp), Hl.reshape(X.shape[0], ncomp))



class TestPCVPCAMethods(unittest.TestCase):

    def setUp(self):
        self.X = np.array([1, 3, 17, 19, 10, 14, 7, 13, 9, 11, 2, 15, 8, 6, 4, 12,
         22, 18, 4, 14, 13, 12, 1, 15]).reshape(3, 8).transpose()



    def test_case1(self):
        """
        Test first case - all possible components + full cross-validation.
        """

        ncomp = 3
        cv = {'type': 'loo'}
        center = True
        scale = False

        test_pcacase(self.X, ncomp, cv, center, scale,
            T = [
                -10.0140380, -8.5235218, -3.2171212,
                -8.6658232, -4.0965210, 1.7038996,
                11.5878439, -1.8647639, -3.9084240,
                0.7731364, 9.8942609, 4.5041217,
                -0.6373488, -0.3702545, -0.4877993,
                2.8319600, 1.4993523, -2.7862511,
                5.6623197, -9.4108933, 6.3367557,
                -1.7167852, 4.7241645, 1.0080517
            ],
            H = [
                1.46251181, 3.58303765,  5.59816235,
                1.09521796, 1.58503715,  2.15030603,
                1.95833139, 2.05982837,  5.03403003,
                0.00871752, 2.86612471,  6.81603569,
                0.00592428, 0.00992562,  0.05625425,
                0.11696479, 0.18258125,  1.69407924,
                0.46759481, 3.05263393, 10.87073001,
                0.04298460, 0.69439613,  0.89224461
            ],
            Q = [
                83.0002919, 10.3498690,  0,
                19.6847583,  2.9032738,  0,
                18.7531225, 15.2757780,  0,
                118.1835101, 20.2871123,  0,
                0.3750365,  0.2379482,  0,
                10.0112525,  7.7631951,  0,
                128.7193857, 40.1544729,  0,
                23.3338984,  1.0161683,  0
            ])



    def test_case2(self):
        """
        Test second case - 2 components + ven with 3 segments.
        """

        ncomp = 2
        cv = {'type': 'ven', 'nseg': 3}
        center = True
        scale = False

        test_pcacase(self.X, ncomp, cv, center, scale,
            T = [
                -12.2937342,  -3.5829202,
                -8.8832034,  -3.6204597,
                11.5446749,  -1.5220728,
                1.1951282, 10.6782087,
                -0.6186727,  -0.4293741,
                2.6220429,  1.3324750,
                7.0616804,  -8.9658336,
                -2.2733024,  4.4926312
            ],
            H = [
                2.20418755, 2.57888410,
                1.15085371, 1.53344303,
                1.94376754, 2.01138772,
                0.02083099, 3.34897591,
                0.00558217, 0.01096334,
                0.10026757, 0.15209070,
                0.72727227, 3.07358973,
                0.07536945, 0.66449378
            ],
            Q = [
                32.1453500, 19.308033,
                15.8699472,  2.762219,
                19.7517308, 17.435025,
                117.3529185,  3.328778,
                0.3984941,  0.214132,
                11.1561408,  9.380651,
                110.9139196, 30.527747,
                21.1133461,  0.929611
            ])



    def test_case3(self):
        """
        Test third case - 1 component + ven with 2 segments.
        """

        ncomp = 1
        cv = {'type': 'ven', 'nseg': 2}
        center = True
        scale = False

        test_pcacase(self.X, ncomp, cv, center, scale,
            T = [
                -9.39633799,
                -8.75468266,
                4.50645582,
                0.01524855,
                -0.64466043,
                2.54815594,
                -3.92591475,
                -2.36235379
            ],
            H = [
                1.28765100,
                1.11779387,
                0.29617687,
                0.00000339,
                0.00606098,
                0.09469628,
                0.22478257,
                0.08138995
            ],
            Q = [
                94.9900823,
                18.1367816,
                132.7231059,
                118.7810175,
                0.3656629,
                11.5381513,
                145.3684434,
                20.7005346
        ])



    def test_reference_cases(self):
        """
        Run several tests by combining PCV parameters and compare outcomes with
        reference data from R.
        """

        # check if reference files exist
        if  not os.path.exists('../.tests') or \
            not os.path.exists('../.tests/data') or \
            not os.path.exists('../.tests/data/corn.csv') or \
            not os.path.exists('../.tests/pcvpca'):

            print('can not find reference files, skipping.')
            return

        X = np.genfromtxt('../.tests/data/corn.csv', delimiter=',')[:, 1:]

        cv_cases = [{'type':'loo'}, {'type':'ven', 'nseg': 4}, {'type':'ven', 'nseg': 10}]
        ncomp_cases = [1, 10, 20, 30]
        scale_cases = [True, False]

        # loop over all combinations of the parameters
        all_cases = list(itertools.product(ncomp_cases, cv_cases, scale_cases))
        for ncomp, cv, scale in all_cases:
            test_pcacase_ref(X, ncomp, cv = cv, center = True, scale = scale)



class TestXpvOrthMethods(unittest.TestCase):

    def test_getxpvorth(self):
        """
        Test of method which compute orthogonal part of Xpv.
        """
        Xk = np.array([[1., 2], [3, 4], [5, 6], [7, 8]])
        qk = np.array([1., 2, 4, 6])
        P = np.array([math.sqrt(2), -math.sqrt(2)])
        RPM = np.eye(2) - np.dot(P, np.transpose(P))
        Y = get_xpvorth(qk, Xk, RPM)

        self.assertEqual(Y.shape, Xk.shape)
        np.testing.assert_array_almost_equal((Y * Y).sum(axis = 1), qk)


if __name__ == '__main__':
    unittest.main()