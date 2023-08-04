import numpy as np
import unittest
import math

from src.prcv.misc import get_cvsettings
from src.prcv.methods import pcvpcr, get_pcamodel
from .common import test_reference_cases


################################
# Helper functions             #
################################

def pcr_predict(X: np.ndarray, Y: np.ndarray, Xpv: np.ndarray, ncomp: int, center: bool = True, scale: bool = False):
    """
    Create a global PCR model using dataset {X, Y}, applies this model to dataset "Xpv"
    and return matrices with predicted Y-values (Yp) and prediction errors (RMSE).
    """

    (nrows, ncols) = X.shape

    # correct shape of Y just in case
    if Y.shape is (nrows, ):
        Y.reshape(nrows, 1)

    nresp = 1

    # compute global mean and standard deviation and autoscale the whole data
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    mY = Y.mean(axis = 0) if center else np.zeros(nresp)
    sY = Y.std(axis = 0, ddof = 1) if scale else np.ones(nresp)

    X = (X - mX) / sX
    Y = (Y - mY) / sY

    # create a global model
    P, s = get_pcamodel(X, ncomp)
    T = np.dot(X, P)
    C = np.dot(np.transpose(Y), T) / (s * s)

    # apply the model
    Xpv = (Xpv - mX) / sX
    Tpv = np.dot(Xpv, P)
    Upv = Tpv / (s / math.sqrt(nrows - 1))

    Ypv = np.zeros((nrows, ncomp))
    Hpv = np.zeros((nrows, ncomp))
    Qpv = np.zeros((nrows, ncomp))

    for a in range(1, ncomp + 1):
        Pa = P[..., :a]
        Ca = C[..., :a]
        Tpva = Tpv[..., :a]
        Upva = Upv[..., :a]
        Epva = Xpv - np.dot(Tpva, np.transpose(Pa))

        Qpv[..., a - 1] = (Epva * Epva).sum(axis = 1)
        Hpv[..., a - 1] = (Upva * Upva).sum(axis = 1)
        Ypv[..., a - 1] = np.dot(Tpva, np.transpose(Ca)).flatten()

    Ypv = Ypv * sY + mY
    Epv = Ypv - Y

    return {'Yp': Ypv, 'T': Tpv, 'H': Hpv, 'Q': Qpv, 'RMSE': np.sqrt((Epv * Epv).mean(axis = 0))}


def test_pcrcase(X, Y, ncomp, cv, center, scale, D, T, Yp):
    """
    Create PCR global model, apply it to generated Xpv set and test equality of
    main outcomes (predictions and errors) entered manually from R tests.
    """
    cvind, cvncomp, cvnseg = get_cvsettings(cv, X.shape[0], ncomp, Y[..., :1])
    Xpv, Dpv = pcvpcr(X, Y, ncomp, center = center, scale = scale, cv = cv)
    r = pcr_predict(X, Y, Xpv, cvncomp, center = center, scale = scale)
    np.testing.assert_array_almost_equal(Dpv, np.array(D).reshape(cvnseg, cvncomp), decimal = 5)
    np.testing.assert_array_almost_equal(r['Yp'], np.array(Yp).reshape(X.shape[0], cvncomp), decimal = 5)
    np.testing.assert_array_almost_equal(r['T'], np.array(T).reshape(X.shape[0], cvncomp), decimal = 5)


def test_pcrcase_ref(X, Y, ncomp, cv, center, scale, path, file_suffix):
    """
    Create PCR global model, apply it to generated Xpv set and test equality of
    main outcomes (predictions and scalars) taken from a
    reference file generated from R tests
    """

    # read reference values for predicted responses (global and local scope)
    Ypvg = np.genfromtxt(path + 'Ypvg' + file_suffix, delimiter=',', ndmin = 2)
    Ypvl = np.genfromtxt(path + 'Ypvl' + file_suffix, delimiter=',', ndmin = 2)

    # read reference values for matrix with scalars (global and local scope)
    Dg = np.genfromtxt(path + 'Dg' + file_suffix, delimiter=',', ndmin = 2)
    Dl = np.genfromtxt(path + 'Dl' + file_suffix, delimiter=',', ndmin = 2)

    # compute PV-sets for global and local scope
    Xpvg, Dpvg = pcvpcr(X, Y, ncomp, center = center, scale = scale, cv = cv, cvscope = 'global')
    Xpvl, Dpvl = pcvpcr(X, Y, ncomp, center = center, scale = scale, cv = cv, cvscope = 'local')

    # make predictions for each PV-set
    rg = pcr_predict(X, Y, Xpvg, ncomp, center = center, scale = scale)
    rl = pcr_predict(X, Y, Xpvl, ncomp, center = center, scale = scale)

    # compare scalars with reference values
    np.testing.assert_array_almost_equal(Dpvg, Dg)
    np.testing.assert_array_almost_equal(Dpvl, Dl)

    # compare predictions with reference values
    np.testing.assert_array_almost_equal(rg['Yp'], Ypvg)
    np.testing.assert_array_almost_equal(rl['Yp'], Ypvl)



################################
# Tests                        #
################################

class TestPCVPCRMethods(unittest.TestCase):

    def setUp(self):
        """ Set up simple data for manual test. """
        self.Y = np.array([32, 35, 36, 37, 42, 43, 43, 44]).reshape(8, 1)
        self.X = np.array([
            150, 41, 28000, 119,
            160, 48, 31000, 129,
            166, 47, 28000, 112,
            166, 49, 14000, 123,
            175, 67, 38000, 105,
            180, 80, 30000, 129,
            181, 75, 31000, 105,
            180, 81, 42000, 113
        ]).reshape(8, 4)


    def test_manual(self):
        """ Test with manually entered reference values from R. """

        ncomp = 3
        cv = {'type': 'ven', 'nseg': 4}
        center = True
        scale = True

        test_pcrcase(self.X, self.Y, ncomp = ncomp, cv = cv, center=center, scale=scale,
            D = [
                0.8928038, 0.6953515,  1.7268124,
                0.9272071, 1.4768482, -0.4511441,
                0.9676971, 1.1075906,  1.1945457,
                1.1402709, 0.6762249, -0.7272937
            ],
            T = [
                -1.5959914, -0.35669193,  -1.7114956,
                -1.1768841, -0.90612631, 0.4232573,
                -0.6259745, -0.77332366, 0.4885416,
                -1.3751630,  0.30118095, -1.3175077,
                1.2243539, -0.61895968,  -0.4291630,
                0.5703641,  1.67461147, 0.6487102,
                1.4017631, -0.12549636, 1.2603507,
                1.9661328,  0.07996545, 0.8937680
            ],
            Yp = [
                34.74990, 34.24415, 32.50785,
                35.86598, 34.58119, 35.01058,
                37.33304, 36.23655, 36.73217,
                35.33797, 35.76501, 34.42841,
                42.26043, 41.38281, 40.94743,
                40.51887, 42.89329, 43.55140,
                42.73287, 42.55493, 43.83354,
                44.23578, 44.34916, 45.25588
            ]
        )


    def test_references(self):
        """ Test combination of settings by comparing with reference values. """
        test_reference_cases('pcvpcr', test_pcrcase_ref)

if __name__ == '__main__':
    unittest.main()