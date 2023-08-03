import numpy as np

import unittest
import itertools
import math
import os

from src.pcv.misc import get_cvsettings
from src.pcv.methods import pcvpls, simpls

def scale(X:np.ndarray) -> np.ndarray:
    return (X - X.mean(axis = 0)) / X.std(axis = 0, ddof = 1)

def pls_predict(X: np.ndarray, Y: np.ndarray, Xpv: np.ndarray, ncomp: int, center: bool = True, scale: bool = False):
    """
    Create a global PLS model using dataset {X, Y}, applies this model to dataset "Xpv"
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
    R, P, C = simpls(X, Y, ncomp)

    # apply the model
    Xpv = (Xpv - mX) / sX
    Tpv = np.dot(Xpv, R)

    Ypv = np.zeros((nrows, ncomp))
    for a in range(1, ncomp + 1):
        Ca = C[..., :a]
        Tpva = Tpv[..., :a]
        Ypv[..., a - 1] = np.dot(Tpva, np.transpose(Ca)).flatten()

    Ypv = Ypv * sY + mY
    Epv = Ypv - Y

    return {'Yp': Ypv, 'T': Tpv, 'RMSE': np.sqrt((Epv * Epv).mean(axis = 0))}



def test_plscase(X, Y, ncomp, cv, center, scale, D, T, Yp):
    """
    Create PLS global model, apply it to generated Xpv set and test equality of
    main outcomes (predictions and errors) entered manually from R tests.
    """
    cvind, cvncomp, cvnseg = get_cvsettings(cv, X.shape[0], ncomp, Y[..., :1])
    Xpv, Dpv = pcvpls(X, Y, ncomp, center = center, scale = scale, cv = cv)
    r = pls_predict(X, Y, Xpv, cvncomp, center = center, scale = scale)
    np.testing.assert_array_almost_equal(Dpv, np.array(D).reshape(cvnseg, cvncomp), decimal = 5)
    np.testing.assert_array_almost_equal(r['Yp'], np.array(Yp).reshape(X.shape[0], cvncomp), decimal = 5)
    np.testing.assert_array_almost_equal(r['T'], np.array(T).reshape(X.shape[0], cvncomp), decimal = 5)



def test_plscase_ref(X, Y, ncomp, cv, center, scale):
    """
    Create PLS global model, apply it to generated Xpv set and test equality of
    main outcomes (predictions and scalars) taken from a
    reference file generated from R tests
    """

    file_suffix = '-' + str(ncomp) + '-' + str(scale).upper() + '-' + \
        (cv['type'] + str(cv['nseg']) if 'nseg' in cv else cv['type']) + '.csv'

    Ypvg = np.genfromtxt('../.tests/pcvpls/Ypvg' + file_suffix, delimiter=',')
    Dg = np.genfromtxt('../.tests/pcvpls/Dg' + file_suffix, delimiter=',')
    Ypvl = np.genfromtxt('../.tests/pcvpls/Ypvl' + file_suffix, delimiter=',')
    Dl = np.genfromtxt('../.tests/pcvpls/Dl' + file_suffix, delimiter=',')

    Xpvg, Dpvg = pcvpls(X, Y, ncomp, center = center, scale = scale, cv = cv, cvscope = 'global')
    Xpvl, Dpvl = pcvpls(X, Y, ncomp, center = center, scale = scale, cv = cv, cvscope = 'local')

    rg = pls_predict(X, Y, Xpvg, ncomp, center = center, scale = scale)
    rl = pls_predict(X, Y, Xpvl, ncomp, center = center, scale = scale)

    np.testing.assert_array_almost_equal(Dpvg.flatten(), Dg.flatten())
    np.testing.assert_array_almost_equal(Dpvl.flatten(), Dl.flatten())
    np.testing.assert_array_almost_equal(rg['Yp'].flatten(), Ypvg.flatten())
    np.testing.assert_array_almost_equal(rl['Yp'].flatten(), Ypvl.flatten())



class TestPCVPLSMethods(unittest.TestCase):

    def setUp(self):
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


    def test_simpls(self):
        """ Tests for SIMPLS implementation """
        R, P, C = simpls(scale(self.X), scale(self.Y), 3)
        np.testing.assert_array_almost_equal(C, np.array([
            [2.559598, 0.5997323, 0.1123446]
        ]))
        np.testing.assert_array_almost_equal(R, np.array([
            [0.15492668, 0.1881235, 0.12558132],
            [0.15317320, 0.1316393, -0.06823508],
            [0.08218702, -0.2824183, -0.32337181],
            [-0.05900391, 0.2131759, -0.34074718]
        ]))
        np.testing.assert_array_almost_equal(P, np.array([
            [ 2.481456,  0.7967174, 0.3398910],
            [ 2.526248,  0.6219986, -0.3463264],
            [ 1.828851, -1.5105447, -1.1672083],
            [-1.326946,  1.6025908, -1.6324190]
        ]))


    def test_manual(self):
        """
        Test with manually entered reference values from R.
        """

        ncomp = 3
        cv = {'type': 'ven', 'nseg': 4}
        center = True
        scale = True

        test_plscase(self.X, self.Y, ncomp = ncomp, cv = cv, center=center, scale=scale,
            D = [
                0.7618305, 0.9701346,  0.8583589,
                0.8664892, 1.1233015, -0.4189311,
                0.9040871, 1.0111721,  1.4835077,
                0.9153578, 0.4453001, -1.0219994,
            ],
            T = [
                -0.4395929, -0.47753680, -0.1098689,
                -0.2967576, -0.25378922,  0.2791594,
                -0.1757461, -0.17091191,  0.5104952,
                -0.3086492,  0.03666759, -2.1360181,
                0.2527472, -0.40208568,  0.2000121,
                0.1879328,  0.80050607,  0.3795580,
                0.3488912,  0.11874715,  0.8917092,
                0.4663885,  0.04006119,  1.4380749,
            ],
            Yp = [
                33.89666, 32.59770, 32.54171,
                35.55487, 34.86453, 35.00677,
                36.95972, 36.49482, 36.75494,
                35.41682, 35.51656, 34.42815,
                41.93420, 40.84048, 40.94240,
                41.18176, 43.35924, 43.55264,
                43.05036, 43.37337, 43.82774,
                44.41442, 44.52339, 45.25616,
            ]
        )



    def test_reference_cases(self):
        """
        Run several tests by combining PCV parameters and compare outcomes with
        reference data from R.
        """

        # check if reference files exist
        if  not os.path.exists('../.tests') or \
            not os.path.exists('../.tests/data') or \
            not os.path.exists('../.tests/data/corn.csv') or \
            not os.path.exists('../.tests/pcvpls'):

            print('can not find reference files, skipping.')
            return


        D = np.genfromtxt('../.tests/data/corn.csv', delimiter=',')
        X = D[:, 1:]
        Y = D[:, :1]

        cv_cases = [{'type':'loo'}, {'type':'ven', 'nseg': 4}, {'type':'ven', 'nseg': 10}]
        ncomp_cases = [1, 10, 20, 30]
        scale_cases = [True, False]

        # loop over all combinations of the parameters
        all_cases = list(itertools.product(ncomp_cases, cv_cases, scale_cases))
        for ncomp, cv, scale in all_cases:
            test_plscase_ref(X, Y, ncomp, cv = cv, center = True, scale = scale)



if __name__ == '__main__':
    unittest.main()