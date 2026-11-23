import itertools
import os
import unittest
import warnings

import numpy as np

from src.prcv.methods import pcvpca
from src.prcv.misc import get_cvsettings
from src.prcv.pcv3w import pcv3wfu, pcv3wbu, pcv3whsv, _align3w, _complementdir3w, \
    _getxpvorth3w, _proj3w



################################
# Helper functions             #
################################

def simdata3w(I: int = 24, J: int = 10, K: int = 8, seed: int = 1) -> np.ndarray:
    """
    Simulated three-way data: three trilinear components plus noise.
    """
    rng = np.random.default_rng(seed)
    A = rng.normal(size = (I, 3))
    B = rng.uniform(size = (J, 3))
    C = rng.uniform(size = (K, 3))
    X = np.einsum('if,jf,kf->ijk', A, B, C)
    return X + rng.normal(scale = 0.05, size = X.shape) + 5



def ref_prep(X: np.ndarray, center: bool, scale: str):
    """
    Reference preprocessing, written independently of the package code.
    Returns preprocessed X and a function which applies the same preprocessing to other data.
    """
    (_, J, K) = X.shape
    M = X.mean(axis = 0) if center else np.zeros((J, K))
    Xp = X - M
    S2 = np.array([Xp[:, j, :].std(ddof = 1) for j in range(J)]) if scale in ('mode2', 'both') else np.ones(J)
    Xp = Xp / S2[None, :, None]
    S3 = np.array([Xp[:, :, k].std(ddof = 1) for k in range(K)]) if scale in ('mode3', 'both') else np.ones(K)
    Xp = Xp / S3[None, None, :]
    return (Xp, lambda Y: (Y - M) / S2[None, :, None] / S3[None, None, :])



def ref_factors(X: np.ndarray, mode: int, R: int) -> np.ndarray:
    """
    Left singular vectors of mode 2 or mode 3 unfolding.
    """
    M = np.concatenate([X[i] for i in range(X.shape[0])], axis = 1) if mode == 2 else \
        np.concatenate([X[i].T for i in range(X.shape[0])], axis = 1)
    return np.linalg.svd(M, full_matrices = False)[0][:, :R]



def is_optimal(M: np.ndarray, perm: np.ndarray) -> bool:
    """
    True if assignment perm (numbered from 1) gives the maximal sum of M[i, perm[i]]
    (True without checking above 6 x 6).
    """
    n = M.shape[0]
    if n > 6:
        return True
    best = max(M[np.arange(n), list(p)].sum() for p in itertools.permutations(range(n)))
    return abs(M[np.arange(n), perm - 1].sum() - best) < 1e-12



def check_hsv_cv(test: unittest.TestCase, X: np.ndarray, Xpv: np.ndarray, info: dict,
        center: bool, scale: str, tol: float = 1e-8):
    """
    Check that projections of HSV PV-set on the global model equal cross-validation results.
    """
    (R2, R3) = info['ncomp']
    cvind = info['cvind']
    al = info['alignment']

    (Xp, apply) = ref_prep(X, center, scale)
    Pp = apply(Xpv)
    U2 = ref_factors(Xp, 2, R2)
    U3 = ref_factors(Xp, 3, R3)

    for k in np.unique(cvind):
        L2 = ref_factors(Xp[cvind != k], 2, R2)
        L3 = ref_factors(Xp[cvind != k], 3, R3)

        # reported matching must reach the maximal sum of absolute cosines (checked up to 6 x 6)
        p2 = al[(al[:, 0] == k) & (al[:, 1] == 2), 3].astype(int)
        p3 = al[(al[:, 0] == k) & (al[:, 1] == 3), 3].astype(int)
        if info['align'] == 'optimal':
            test.assertTrue(is_optimal(np.abs(U2.T @ L2), p2))
            test.assertTrue(is_optimal(np.abs(U3.T @ L3), p3))

        # local factors in the order reported by the method, signs matched to the global ones
        U2k = L2[:, p2 - 1]
        U3k = L3[:, p3 - 1]
        U2k = U2k * np.sign((U2 * U2k).sum(axis = 0))
        U3k = U3k * np.sign((U3 * U3k).sum(axis = 0))

        for i in np.flatnonzero(cvind == k):
            Xi = Xp[i]
            Pi = Pp[i]

            # all components: scores and orthogonal distance
            Tcv = U2k.T @ Xi @ U3k
            Tpv = U2.T @ Pi @ U3
            np.testing.assert_allclose(Tpv, Tcv, rtol = tol, atol = tol * np.abs(Tcv).max())
            qcv = ((Xi - U2k @ Tcv @ U3k.T) ** 2).sum()
            qpv = (Pi ** 2).sum() - (Tpv ** 2).sum()
            test.assertAlmostEqual(qpv, qcv, delta = tol * (Xi ** 2).sum())

            # valid smaller numbers of components: compare with local models fitted directly
            for r2 in info['valid_subranks']['r2']:
                for r3 in info['valid_subranks']['r3']:
                    l2 = L2[:, :r2]
                    l3 = L3[:, :r3]
                    qcv = ((Xi - l2 @ (l2.T @ Xi @ l3) @ l3.T) ** 2).sum()
                    qpv = (Pi ** 2).sum() - (Tpv[:r2, :r3] ** 2).sum()
                    test.assertAlmostEqual(qpv, qcv, delta = tol * (Xi ** 2).sum())



def ref_qh(X: np.ndarray, Xpv: np.ndarray, R: tuple, center: bool, scale: str):
    """
    q and h of a PV-set on the global model (h uses the pseudo-inverse of the score scatter).
    """
    (Xp, apply) = ref_prep(X, center, scale)
    U2 = ref_factors(Xp, 2, R[0])
    U3 = ref_factors(Xp, 3, R[1])

    # scores are vectorized column by column, as in R
    scores = lambda A: np.array([(U2.T @ A[i] @ U3).flatten(order = 'F') for i in range(A.shape[0])])
    _, s, Vt = np.linalg.svd(scores(Xp), full_matrices = False)
    use = s > np.sqrt(np.finfo(float).eps) * s[0]

    Pp = apply(Xpv)
    T = scores(Pp)
    q = (Pp ** 2).sum(axis = (1, 2)) - (T ** 2).sum(axis = 1)
    h = (((T @ Vt[use].T) / s[use]) ** 2).sum(axis = 1)
    return (q, h)



################################
# Tests                        #
################################

class TestPCV3WHSV(unittest.TestCase):

    def test_preprocessing(self):
        """ HSV PV-set reproduces cross-validation for all preprocessing settings. """
        X = simdata3w()
        for center in (True, False):
            for scale in ('none', 'mode2', 'mode3', 'both'):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    Xpv, info = pcv3whsv(X, (3, 3), center = center, scale = scale, seed = 1)
                self.assertEqual(Xpv.shape, X.shape)
                self.assertEqual(info['ncomp'], (3, 3))
                check_hsv_cv(self, X, Xpv, info, center, scale)

    def test_splits_and_options(self):
        """ HSV PV-set reproduces cross-validation for different splits and settings. """
        X = simdata3w(seed = 2)
        cvs = [{'type': 'ven', 'nseg': 3}, {'type': 'rand', 'nseg': 5}, {'type': 'loo'},
            np.tile(np.arange(1, 5), 6)]
        for cv in cvs:
            for orth in ('random', 'canonical'):
                for align in ('optimal', 'greedy'):
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        Xpv, info = pcv3whsv(X, (4, 3), cv = cv, orth = orth, align = align, seed = 3)
                    check_hsv_cv(self, X, Xpv, info, True, 'none')

    def test_subranks(self):
        """ Some smaller numbers of components are not valid: warning, and they really fail. """
        X = simdata3w(seed = 9)
        with warnings.catch_warnings(record = True) as w:
            warnings.simplefilter('always')
            Xpv, info = pcv3whsv(X, (5, 5), seed = 1)
        vs = info['valid_subranks']
        self.assertIn(5, vs['r2'])
        self.assertIn(5, vs['r3'])
        self.assertLess(len(vs['r2']), 5)
        self.assertTrue(any('smaller numbers of components' in str(x.message) for x in w))
        check_hsv_cv(self, X, Xpv, info, True, 'none')

        # omitted numbers of components really do not reproduce cross-validation
        cvind = info['cvind']
        M = X.mean(axis = 0)
        Xc = X - M
        Pc = Xpv - M
        U2 = ref_factors(Xc, 2, 5)
        U3 = ref_factors(Xc, 3, 5)
        for mode in (2, 3):
            for r in sorted(set(range(1, 6)) - set(vs['r2' if mode == 2 else 'r3'])):
                r2 = r if mode == 2 else 5
                r3 = r if mode == 3 else 5
                maxdiff = 0
                for k in np.unique(cvind):
                    l2 = ref_factors(Xc[cvind != k], 2, r2)
                    l3 = ref_factors(Xc[cvind != k], 3, r3)
                    for i in np.flatnonzero(cvind == k):
                        qcv = ((Xc[i] - l2 @ (l2.T @ Xc[i] @ l3) @ l3.T) ** 2).sum()
                        qpv = (Pc[i] ** 2).sum() - ((U2[:, :r2].T @ Pc[i] @ U3[:, :r3]) ** 2).sum()
                        maxdiff = max(maxdiff, abs(qpv - qcv) / qcv)
                self.assertGreater(maxdiff, 1e-6, f'mode {mode}, {r} components')

    def test_info(self):
        """ Information returned together with the PV-set is correct. """
        X = simdata3w()
        Xpv, info = pcv3whsv(X, (3, 2), scale = True, orth = 'canonical', seed = 1)
        np.testing.assert_array_equal(info['cvind'], np.tile(np.arange(1, 5), 6))
        self.assertEqual(info['orth'], 'canonical')
        self.assertEqual(info['align'], 'optimal')
        self.assertEqual(info['scale'], 'both')
        self.assertEqual(info['prep']['M'].shape, (10, 8))
        self.assertEqual(info['alignment'].shape, (4 * (3 + 2), 5))
        self.assertTrue(np.all(info['alignment'][:, 4] <= 1 + 1e-12))

    def test_canonical_is_deterministic(self):
        """ Canonical orthogonal part does not depend on random numbers. """
        X = simdata3w()
        Xpv1, _ = pcv3whsv(X, (3, 3), orth = 'canonical', seed = 1)
        Xpv2, _ = pcv3whsv(X, (3, 3), orth = 'canonical', seed = 2)
        np.testing.assert_array_equal(Xpv1, Xpv2)

    def test_ncomp_reduced(self):
        """ Too large number of components is reduced with a warning. """
        X = simdata3w(I = 12, J = 4, K = 5)
        with self.assertWarnsRegex(UserWarning, 'reduced'):
            Xpv, info = pcv3whsv(X, (6, 6), cv = {'type': 'ven', 'nseg': 3}, seed = 1)
        self.assertEqual(info['ncomp'], (4, 5))
        check_hsv_cv(self, X, Xpv, info, True, 'none')

    def test_alignment(self):
        """ Optimal and greedy matching give expected assignments. """
        M = np.array([[0.90, 0.80, 0.10], [0.85, 0.10, 0.30], [0.20, 0.75, 0.70]])
        Ug = np.eye(3)
        Ul = M * np.array([1, -1, 1])

        Ua, perm, cos = _align3w(Ug, Ul, 'optimal')
        np.testing.assert_array_equal(perm, [2, 1, 3])
        np.testing.assert_allclose(cos, [0.80, 0.85, 0.70])
        self.assertTrue(np.all((Ug * Ua).sum(axis = 0) > 0))

        _, perm, _ = _align3w(Ug, Ul, 'greedy')
        np.testing.assert_array_equal(perm, [1, 3, 2])

    def test_fixed_direction(self):
        """ Orthogonal part uses the fixed direction when needed. """
        rng = np.random.default_rng(5)
        U2 = np.linalg.qr(rng.normal(size = (6, 2)))[0]
        U3 = np.linalg.qr(rng.normal(size = (5, 2)))[0]

        # direction does not depend on signs of the factors
        cdir = _complementdir3w(U2, U3)
        np.testing.assert_allclose(_complementdir3w(-U2, -U3), cdir)
        self.assertAlmostEqual((cdir ** 2).sum(), 1.0)
        self.assertAlmostEqual((_proj3w(cdir, U2, U3) ** 2).sum(), 0.0)

        # residual lying inside the model space: the fixed direction is used with the right norm
        E = U2 @ np.array([[1.0, 3.0], [2.0, 4.0]]) @ U3.T
        for orth in ('canonical', 'random'):
            Xk = [E] if orth == 'canonical' else [E, 2 * E]
            Eres = [E] if orth == 'canonical' else [E, E]
            Eo = _getxpvorth3w(Xk, Eres, U2, U3, cdir, orth)
            np.testing.assert_allclose(Eo[0], cdir * np.sqrt((E ** 2).sum()))

        # full model space in both modes: no direction
        self.assertIsNone(_complementdir3w(np.eye(3), np.eye(2)))



class TestPCV3WUnfolding(unittest.TestCase):

    def test_fu(self):
        """ FU PV-set is the same as pcvpca() on unfolded data. """
        X = simdata3w()
        (I, J, K) = X.shape
        for scale in (False, True):
            Xpv = pcv3wfu(X, ncomp = 6, scale = scale, seed = 5)
            np.random.seed(5)
            Xpv2 = pcvpca(X.reshape(I, J * K), 6, scale = scale)
            np.testing.assert_array_equal(Xpv, Xpv2.reshape(I, J, K))

            # column (j, k) of the unfolded data is column j * K + k
            for j in (0, J - 1):
                for k in (0, 1, K - 1):
                    np.testing.assert_array_equal(Xpv[:, j, k], Xpv2[:, j * K + k])

    def test_bu(self):
        """ BU PV-set is the same as pcvpca() on each block with one split. """
        X = simdata3w()
        (I, J, K) = X.shape
        cv = {'type': 'rand', 'nseg': 4}
        for mode in (2, 3):
            Xpv = pcv3wbu(X, mode = mode, ncomp = 4, cv = cv, seed = 5)
            np.random.seed(5)
            cvind = get_cvsettings(cv, I, 1)[0]
            for b in range(J if mode == 2 else K):
                Xb = X[:, b, :] if mode == 2 else X[:, :, b]
                Pb = Xpv[:, b, :] if mode == 2 else Xpv[:, :, b]
                np.testing.assert_array_equal(Pb, pcvpca(Xb, 4, cv = cvind))



class TestPCV3WSeedAndInput(unittest.TestCase):

    def test_seed(self):
        """ Seed makes results reproducible and keeps the caller's random state. """
        X = simdata3w()
        funs = [
            lambda s: pcv3whsv(X, (3, 3), seed = s)[0],
            lambda s: pcv3wfu(X, 5, seed = s),
            lambda s: pcv3wbu(X, 2, 4, seed = s)
        ]
        for f in funs:
            np.random.seed(99)
            r1 = np.random.uniform()
            np.random.seed(99)
            P1 = f(10)
            r2 = np.random.uniform()
            self.assertEqual(r1, r2)
            np.testing.assert_array_equal(P1, f(10))
            self.assertFalse(np.allclose(P1, f(11)))

    def test_wrong_input(self):
        """ Wrong input gives an error. """
        X = simdata3w()
        with self.assertRaises(ValueError):
            pcv3whsv(np.ones((4, 5)), (2, 2))
        with self.assertRaises(ValueError):
            pcv3whsv(X[:, :, :1], (2, 2))
        for ncomp in (2, (0, 2), (2.5, 2), (2, 2, 2), 'ab', None):
            with self.assertRaisesRegex(ValueError, 'two positive integers'):
                pcv3whsv(X, ncomp)

        Xna = X.copy()
        Xna[0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            pcv3whsv(Xna, (2, 2))

        Xc = X.copy()
        Xc[:, 1, :] = 1.0
        with self.assertRaisesRegex(ValueError, 'zero variance'):
            pcv3whsv(Xc, (2, 2), scale = 'mode2')

        with self.assertRaises(ValueError):
            pcv3whsv(X, (2, 2), scale = 'wrong')
        with self.assertRaises(ValueError):
            pcv3whsv(X, (2, 2), orth = 'wrong')
        with self.assertRaises(ValueError):
            pcv3whsv(X, (2, 2), align = 'wrong')
        with self.assertRaises(TypeError):
            pcv3whsv(X, (2, 2), cvscope = 'local')
        with self.assertRaisesRegex(ValueError, 'mode'):
            pcv3wbu(X, mode = 1)
        with self.assertRaisesRegex(ValueError, '3D array'):
            pcv3wfu(np.ones((4, 5)))

    def test_wrong_segments(self):
        """ Wrong segment numbers in cv give an error. """
        X = simdata3w()
        cvs = [np.tile([0, 1, 2], 8), np.tile([1, 3], 12), np.ones(24), np.tile([1.5, 2.0], 12)]
        for cv in cvs:
            with self.assertRaisesRegex(ValueError, 'segment numbers'):
                pcv3whsv(X, (2, 2), cv = cv)
            with self.assertRaisesRegex(ValueError, 'segment numbers'):
                pcv3wbu(X, 2, 3, cv = cv)
            with self.assertRaisesRegex(ValueError, 'segment numbers'):
                pcv3wfu(X, 3, cv = cv)

        # vector of whole numbers given as floats or as a column works
        Xpv1 = pcv3wfu(X, 3, cv = np.tile([1.0, 2.0], 12), seed = 1)
        Xpv2 = pcv3wfu(X, 3, cv = np.tile([1, 2], 12).reshape(-1, 1), seed = 1)
        np.testing.assert_array_equal(Xpv1, Xpv2)

    def test_numpy_bool_scale(self):
        """ NumPy booleans are accepted as value for scale. """
        X = simdata3w()
        _, info = pcv3whsv(X, (2, 2), scale = np.True_, seed = 1)
        self.assertEqual(info['scale'], 'both')



class TestPCV3WReferences(unittest.TestCase):

    def test_references(self):
        """ HSV results match the reference values computed in R. """
        data_file = '../.tests/data/sim3w-20x8x6.csv'
        case_dir = '../.tests/pcv3whsv/'
        if not os.path.exists(data_file) or not os.path.exists(case_dir):
            self.skipTest('reference files not found')

        # dataset is stored unfolded, I x (J*K) with mode 2 running fastest
        X = np.genfromtxt(data_file, delimiter = ',').reshape((20, 8, 6), order = 'F')
        rcsv = lambda name: np.genfromtxt(case_dir + name, delimiter = ',')

        cases = [((3, 3), center, scale, 'optimal', {'type': 'ven', 'nseg': 4})
            for center in (True, False) for scale in ('none', 'mode2', 'mode3', 'both')]
        cases += [((4, 2), True, 'none', align, cv) for align in ('optimal', 'greedy')
            for cv in ({'type': 'ven', 'nseg': 4}, {'type': 'loo'})]

        for (ncomp, center, scale, align, cv) in cases:
            cvtext = cv['type'] + (str(cv['nseg']) if 'nseg' in cv else '')
            sfx = f'-{ncomp[0]}x{ncomp[1]}-{str(center).upper()}-{scale}-{align}-{cvtext}.csv'

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Xpv, info = pcv3whsv(X, ncomp, center = center, scale = scale, cv = cv,
                    orth = 'canonical', align = align)
                Xpvr, infor = pcv3whsv(X, ncomp, center = center, scale = scale, cv = cv,
                    orth = 'random', align = align, seed = 1)

            np.testing.assert_allclose(Xpv.reshape((20, 48), order = 'F'), rcsv('Xpv' + sfx),
                rtol = 1e-6, atol = 1e-8, err_msg = 'Xpv' + sfx)
            np.testing.assert_allclose(info['alignment'], rcsv('align' + sfx), rtol = 1e-6,
                atol = 1e-8, err_msg = 'align' + sfx)

            # random orthogonal part: q and h are still determined (equal to cross-validation)
            (q, h) = ref_qh(X, Xpvr, infor['ncomp'], center, scale)
            np.testing.assert_allclose(q, rcsv('Q' + sfx), rtol = 1e-6, err_msg = 'Q' + sfx)
            np.testing.assert_allclose(h, rcsv('H' + sfx), rtol = 1e-6, err_msg = 'H' + sfx)



class TestPCV3WJuices(unittest.TestCase):

    def test_juices(self):
        """ HSV, FU and BU work on Juices data (private, set PCV_JUICES to the path of juices.csv). """
        path = os.environ.get('PCV_JUICES', '')
        if path == '' or not os.path.exists(path):
            self.skipTest('Juices data not available (PCV_JUICES not set)')

        # first column has sample names, second class names, the rest are values
        D = np.genfromtxt(path, delimiter = ',', skip_header = 1)[:, 2:]
        X = D.reshape((D.shape[0], 81, 61)).transpose(0, 2, 1)

        cv = {'type': 'ven', 'nseg': 3}
        Xpv, info = pcv3whsv(X, (10, 10), cv = cv, seed = 1)
        self.assertEqual(info['ncomp'], (10, 10))
        check_hsv_cv(self, X, Xpv, info, True, 'none', tol = 1e-6)

        Xfu = pcv3wfu(X, 20, cv = cv, seed = 1)
        self.assertEqual(Xfu.shape, X.shape)
        self.assertTrue(np.all(np.isfinite(Xfu)))

        for mode in (2, 3):
            Xbu = pcv3wbu(X, mode, 20, cv = cv, seed = 1)
            self.assertEqual(Xbu.shape, X.shape)
            self.assertTrue(np.all(np.isfinite(Xbu)))



if __name__ == '__main__':
    unittest.main()
