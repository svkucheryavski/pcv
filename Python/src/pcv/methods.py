import numpy as np
from scipy.linalg import svd
import math
from .misc import get_cvsettings

def pcvpls(X: np.ndarray,
        Y: np.ndarray,
        ncomp: int,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        cvscope: str = 'global') -> (np.ndarray, np.ndarray):
    """
    Procrustes cross-validation based on PLS decomposition.

    Parameters:
    -----------
    X : numpy.ndarray
        Predictors for global calibration set as 2D Numpy array.
    Y : numpy.ndarray
        Responses for global calibration set as 2D Numpy array (single column).
    ncomp : int
        Number of components to use for generation (must be larger than expected optimal number).
    center : bool, default = True
        Logical, center or not columns of X and Y (better do it).
    scale : bool, default = False
        Logical, standardize or not columns of X and Y.
    cv : array|dict, default = {'type': 'ven', 'nseg', 4}
        Parameters of cross-validation loop, see details below.
    cvscope: str, default = 'global'
        Defines how local calibration set must be centered and/or scaled.

    Returns:
    --------
    Xpv : numpy.ndarray
        2D array of the same size as X with generated PV-set.
    D : numpy.ndarray
        2D array with ck/c scalars for each segment and each component.

    Notes:
    ------
    The method computes Procrustes validation set as matrix Xpv, based on PLS decomposition of
    calibration set {X, Y} and cross-validation resampling. See description of the method in [1].

    Parameter `cv` defines how to split the rows of the training set. The split is similar
    to cross-validation splits, as PCV is based on cross-validation. This parameter can have
    the following values:

    1. A dictionary with 1 or 2 elements: `{'type': 'loo')` or `{'type': 'rand', 'nseg', 4}`. In
    this case `'type'` defines the way to make the splits, you can select one of the following:
    `'loo'` for leave-one-out, `'rand'` for random splits or `'ven'` for Venetian blinds
    (systematic) splits. The second parameter, `'nseg'`, is a number of segments.

    2. A vector (array) with integer numbers, e.g. `np.array([1, 2, 3, 1, 2, 3], dtype = 'int')`.
    In this case number of values in this vector must be the same as number of rows in the training
    set. The values specify which segment a particular row will belong to. In case of the example
    shown here, it is assumed that you have 6 rows in the calibration set, which will be split
    into 3 segments. The first segment will consist of measurements from rows 1 and 4.

    Parameter `cvscope` influences how the Procrustean rule is met. In case of `'global'` scope,
    the rule will be met strictly - error of predictions for PV-set and the global model will be
    identical to the error from conventional cross-validation. In case of `'local'` scope, every
    local model will have its own center and scaling factor and hence the rule will be almost
    met (the errors will be close but not identical).

    References:
    -----------
    1. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev. Procrustes cross-validation of multivariate
    regression models. Analytica Chimica Acta, 1255 (2022) DOI: https://doi.org/10.1016/j.aca.2023.341096
    """

    def get_globalmodel(X: np.ndarray, Y: np.ndarray, ncomp: int) -> dict:
        """ Create a global PLS model. """
        ncols = X.shape[1]
        R, P, C = simpls(X, Y, ncomp)
        PRM = np.eye(ncols) - np.dot(R, P.T)
        return {'R': R, 'P': P, 'C': C, 'PRM': PRM, 'ncomp': ncomp}


    def get_localmodel(Xc: np.ndarray, Yc: np.ndarray, m: dict) -> dict:
        """ Create a local PLS model. """
        ncols = Xc.shape[1]
        Rk, Pk, Ck = simpls(Xc, Yc, ncomp)
        a = np.arccos((normalize(m['R']) * normalize(Rk)).sum(axis = 0)) < math.pi / 2
        Pk = Pk * (2 * a - 1)
        Rk = Rk * (2 * a - 1)
        Ck = Ck * (2 * a - 1)
        PRM = np.eye(ncols) - np.dot(Rk, Pk.T)
        return {'R': Rk, 'P': Pk, 'C': Ck, 'PRM': PRM}

    def get_xpvhat(m: dict, mk: dict, Xk: np.ndarray):
        """ Compute explained part of Xpv for current segment k. """
        Tk = np.dot(Xk, mk['R'])

        dk = np.zeros((m['ncomp']))
        for a in range(0, m['ncomp']):
            ca = m['C'][..., a]
            cka = mk['C'][..., a]
            dk[a] = np.dot(cka.T, ca) / np.dot(ca.T, ca)

        Tpvk = Tk * dk
        Xpvk_hat = np.dot(Tpvk, m['P'].T)
        return (Xpvk_hat, dk)

    def get_qk(Xk: np.ndarray, mk: dict):
        """ Compute vector with orthogonal distances for local model. """
        Ek = np.dot(Xk, mk['PRM'])
        return (Ek * Ek).sum(axis = 1)

    # define function list for pcvreg() method
    funlist = {
        'get_globalmodel': get_globalmodel,
        'get_localmodel': get_localmodel,
        'get_xpvhat': get_xpvhat,
        'get_qk': get_qk
    }

    return pcvreg(X, Y, ncomp, funlist = funlist, cv = cv, center = center, scale = scale, cvscope = cvscope)


def pcvpcr(X: np.ndarray,
        Y: np.ndarray,
        ncomp: int,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        cvscope: str = 'global') -> (np.ndarray, np.ndarray):
    """
    Procrustes cross-validation based on PCR decomposition.

    Parameters:
    -----------
    X : numpy.ndarray
        Predictors for global calibration set as 2D Numpy array.
    Y : numpy.ndarray
        Responses for global calibration set as 2D Numpy array (single column).
    ncomp : int
        Number of components to use for generation (must be larger than expected optimal number).
    center : bool, default = True
        Logical, center or not columns of X and Y (better do it).
    scale : bool, default = False
        Logical, standardize or not columns of X and Y.
    cv : array|dict, default = {'type': 'ven', 'nseg', 4}
        Parameters of cross-validation loop, see details below.
    cvscope: str, default = 'global'
        Defines how local calibration set must be centered and/or scaled.

    Returns:
    --------
    Xpv : numpy.ndarray
        2D array of the same size as X with generated PV-set.
    D : numpy.ndarray
        2D array with ck/c scalars for each segment and each component.

    Notes:
    ------
    The method computes Procrustes validation set as matrix Xpv, based on PCR decomposition of
    calibration set {X, Y} and cross-validation resampling. See description of the method in [1].

    Parameter `cv` defines how to split the rows of the training set. The split is similar
    to cross-validation splits, as PCV is based on cross-validation. This parameter can have
    the following values:

    1. A dictionary with 1 or 2 elements: `{'type': 'loo')` or `{'type': 'rand', 'nseg', 4}`. In
    this case `'type'` defines the way to make the splits, you can select one of the following:
    `'loo'` for leave-one-out, `'rand'` for random splits or `'ven'` for Venetian blinds
    (systematic) splits. The second parameter, `'nseg'`, is a number of segments.

    2. A vector (array) with integer numbers, e.g. `np.array([1, 2, 3, 1, 2, 3], dtype = 'int')`.
    In this case number of values in this vector must be the same as number of rows in the training
    set. The values specify which segment a particular row will belong to. In case of the example
    shown here, it is assumed that you have 6 rows in the calibration set, which will be split
    into 3 segments. The first segment will consist of measurements from rows 1 and 4.

    Parameter `cvscope` influences how the Procrustean rule is met. In case of `'global'` scope,
    the rule will be met strictly - error of predictions for PV-set and the global model will be
    identical to the error from conventional cross-validation. In case of `'local'` scope, every
    local model will have its own center and scaling factor and hence the rule will be almost
    met (the errors will be close but not identical).

    References:
    -----------
    1. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev. Procrustes cross-validation of multivariate
    regression models. Analytica Chimica Acta, 1255 (2022) DOI: https://doi.org/10.1016/j.aca.2023.341096
    """

    def get_globalmodel(X: np.ndarray, Y: np.ndarray, ncomp: int) -> dict:
        """ Create a global PCR model. """
        ncols = X.shape[1]
        P, s = get_pcamodel(X, ncomp)
        PRM = np.eye(ncols) - np.dot(P, P.T)
        T = np.dot(X, P)
        C = np.dot(Y.T, T) / (s * s)
        return {'P': P, 'C': C, 'PRM': PRM, 'ncomp': ncomp}


    def get_localmodel(Xc: np.ndarray, Yc: np.ndarray, m: dict) -> dict:
        """ Create a local PCR model. """
        ncols = Xc.shape[1]
        Pk, sk = get_pcamodel(Xc, m['ncomp'])
        a = np.arccos((m['P'] * Pk).sum(axis = 0)) < math.pi / 2
        Pk = Pk * (2 * a - 1)
        PRMk = np.eye(ncols) - np.dot(Pk, Pk.T)
        Tc = np.dot(Xc, Pk)
        Ck = np.dot(Yc.T, Tc) / (sk * sk)
        return {'P': Pk, 'C': Ck, 'PRM': PRMk}

    def get_xpvhat(m: dict, mk: dict, Xk: np.ndarray):
        """ Compute explained part of Xpv for current segment k. """
        Tk = np.dot(Xk, mk['P'])
        dk = mk['C'] / m['C']
        Tpvk = Tk * dk
        Xpvk_hat = np.dot(Tpvk, m['P'].T)
        return (Xpvk_hat, dk)

    def get_qk(Xk: np.ndarray, mk: dict):
        """ Compute vector with orthogonal distances for local model. """
        Ek = np.dot(Xk, mk['PRM'])
        return (Ek * Ek).sum(axis = 1)

    # define function list for pcvreg() method
    funlist = {
        'get_globalmodel': get_globalmodel,
        'get_localmodel': get_localmodel,
        'get_xpvhat': get_xpvhat,
        'get_qk': get_qk
    }

    return pcvreg(X, Y, ncomp, funlist = funlist, cv = cv, center = center, scale = scale, cvscope = cvscope)


def pcvreg(X: np.ndarray,
        Y: np.ndarray,
        ncomp: int,
        funlist: dict,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        cvscope: str = 'global') -> (np.ndarray, np.ndarray):
    """
    Procrustes cross-validation for regression methods.

    Parameters:
    -----------
    X : numpy.ndarray
        Predictors for global calibration set as 2D Numpy array.
    Y : numpy.ndarray
        Responses for global calibration set as 2D Numpy array.
    ncomp : int
        Number of components to use for generation (must be larger than expected optimal number).
    funlist :  dict
        Dictionary with main PCV functions for particular implementation (PCV or PLS).
    center : bool, default = True
        Logical, center or not columns of X and Y (better do it).
    scale : bool, default = False
        Logical, standardize or not columns of X and Y.
    cv : array|dict, default = {'type': 'ven', 'nseg': 4}
        Parameters of cross-validation loop, see details below.
    cvscope: str, default = 'global'
        Defines how local calibration set must be centered and/or scaled.

    Returns:
    --------
    Xpv : numpy.ndarray
        2D array of the same size as X with generated PV-set.
    D : numpy.ndarray
        2D array with ck/c scalars for each segment and each component.

    Notes:
    ------
    This is a generic method, use `pcvpls()` or `pcvpcr()` instead.

    """

    (nrows, ncols) = X.shape

    # correct shape of Y just in case
    if Y.shape is (nrows, ):
        Y.reshape(nrows, 1)

    nresp = Y.shape[1]

    # get cross-validation settings
    (cvind, cvncomp, cvnseg) = get_cvsettings(cv, nrows, ncomp, Y[..., :1])

    # compute global mean and standard deviation and autoscale the whole data
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    mY = Y.mean(axis = 0) if center else np.zeros(nresp)
    sY = Y.std(axis = 0, ddof = 1) if scale else np.ones(nresp)

    X = (X - mX) / sX
    Y = (Y - mY) / sY


    # create a global model
    m = funlist['get_globalmodel'](X, Y, cvncomp)

    # prepare empty matrix for pseudo-validation set and matrix of scalars
    Xpv = np.zeros((nrows, ncols))
    D = np.zeros((cvnseg, cvncomp))

    for k in range(1, cvnseg + 1):

        indc = cvind != k
        indk = cvind == k

        Xc = X[indc, ...]
        Yc = Y[indc, ...]
        Xk = X[indk, ...]

        # if cvscope is 'local' autoscale local calibration and validation sets
        if cvscope == 'local':
            mXl = Xc.mean(axis = 0) if center else np.zeros(ncols)
            sXl = Xc.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
            mYl = Yc.mean(axis = 0) if center else np.zeros(nresp)
            sYl = Yc.std(axis = 0, ddof = 1) if scale else np.ones(nresp)
            Xc = (Xc - mXl) / sXl
            Yc = (Yc - mYl) / sYl
            Xk = (Xk - mXl) / sXl

        # get local model
        mk = funlist['get_localmodel'](Xc, Yc, m)

        # compute explained part of PV-set for current segment
        Xpvk_hat, dk = funlist['get_xpvhat'](m, mk, Xk)

        # compute the orthogonal component of Xpv
        Xpvk_orth = get_xpvorth(funlist['get_qk'](Xk, mk), Xk, m['PRM'])

        # create Xpvk as a sum of both parts and save to main array
        Xpv[indk, ...] = Xpvk_hat + Xpvk_orth
        D[k - 1, ...] = dk

    # uscenter and unscale the data using global mean and std and return the result
    return (Xpv * sX + mX, D)


def pcvpca(X: np.ndarray,
        ncomp: int,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        cvscope: str = 'global') -> np.ndarray:
    """
    Procrustes cross-validation based on PCA decomposition.

    Parameters:
    -----------
    X : numpy.ndarray
        Global calibration set as 2D Numpy array.
    ncomp : int
        Number of components to use for generation (must be larger than expected optimal number).
    center : bool, default = True
        Logical, center or not columns of X (better do it).
    scale : bool, default = False
        Logical, standardize or not columns of X.
    cv : array|dict, default = {'type': 'ven', 'nseg', 4}
        Parameters of cross-validation loop, see details below.
    cvscope: str, default = 'global'
        Defines how local calibration set must be centered and/or scaled.

    Returns:
    --------
    Xpv : numpy.ndarray
        2D array of the same size as X with generated PV-set.

    Notes:
    ------
    The method computes Procrustes validation set as matrix Xpv, based on PCA decomposition of
    calibration set X and cross-validation resampling. See description of the method in [1].

    Parameter `cv` defines how to split the rows of the training set. The split is similar
    to cross-validation splits, as PCV is based on cross-validation. This parameter can have
    the following values:

    1. A dictionary with 1 or 2 elements: `{'type': 'loo')` or `{'type': 'rand', 'nseg', 4}`. In
    this case `'type'` defines the way to make the splits, you can select one of the following:
    `'loo'` for leave-one-out, `'rand'` for random splits or `'ven'` for Venetian blinds
    (systematic) splits. The second parameter, `'nseg'`, is a number of segments.

    2. A vector (array) with integer numbers, e.g. `np.array([1, 2, 3, 1, 2, 3], dtype = 'int')`.
    In this case number of values in this vector must be the same as number of rows in the training
    set. The values specify which segment a particular row will belong to. In case of the example
    shown here, it is assumed that you have 6 rows in the calibration set, which will be split
    into 3 segments. The first segment will consist of measurements from rows 1 and 4.

    Parameter `cvscope` influences how the Procrustean rule is met. In case of `'global'` scope,
    the rule will be met strictly - distances for PV-set and the global model will be
    identical to the distances from conventional cross-validation. In case of `'local'` scope, every
    local model will have its own center and scaling factor and hence the rule will be almost
    met (the distances will be close but not identical).

    References:
    -----------
    1. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev. Procrustes cross-validation of multivariate
    regression models. Analytica Chimica Acta, 1255 (2022) DOI: https://doi.org/10.1016/j.aca.2023.341096
    """

    (nrows, ncols) = X.shape

    # get cross-validation settings
    (cvind, cvncomp, cvnseg) = get_cvsettings(cv, nrows, ncomp)

    # compute global mean and standard deviation and autoscale the whole data
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    X = (X - mX) / sX


    # create a global model
    P, s = get_pcamodel(X, cvncomp)
    PRM = np.eye(ncols) - np.dot(P, P.T)

    # prepare empty matrix for pseudo-validation set
    Xpv = np.zeros((nrows, ncols))

    for k in range(1, cvnseg + 1):

        # split data to local calibration and validation sets
        indc = cvind != k
        indk = cvind == k

        Xc = X[indc, ...]
        Xk = X[indk, ...]

        # if cvscope is 'local', autoscale the sets locally
        if cvscope == 'local':
            mXl = Xc.mean(axis = 0) if center else np.zeros(ncols)
            sXl = Xc.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
            Xc = (Xc - mXl) / sXl
            Xk = (Xk - mXl) / sXl

        # get local model and correct directions of the loadings
        Pk, sk = get_pcamodel(Xc, cvncomp)
        a = np.arccos((P * Pk).sum(axis = 0)) < math.pi / 2
        Pk = Pk * (2 * a - 1)

        # get scores and residuals by projection local validation set to the local model
        Tk = np.dot(Xk, Pk)
        Ek = Xk - np.dot(Tk, Pk.T)
        qk = (Ek * Ek).sum(axis = 1)

        # compute the explained part of Xpv
        Xpvk_hat = np.dot(Tk, P.T)

        # compute the residusl part of Xpv
        Xpvk_orth = get_xpvorth(qk, Xk, PRM)

        # create and save the Xpv as a sum of the two parts
        Xpv[indk, ...] = Xpvk_hat + Xpvk_orth

    # uscenter and unscale the data using global mean and std and return the result
    return Xpv * sX + mX


def get_pcamodel(X: np.ndarray, ncomp: int) -> dict:
    """ Create PCA model by using SVD and truncating loadings and singular values. """
    U, s, V = np.linalg.svd(X, full_matrices=False)
    return (V[:ncomp, ...].T, s[:ncomp])


def get_xpvorth(qk: np.ndarray, Xk: np.ndarray, PRM: np.ndarray) -> np.ndarray:
    """
    Generate the orthogonal (residual) part for Xpv.

    Parameters:
    -----------
    qk : numpy.ndarray
        Vector with orthogonal distances for cross-validation set for segment k.
    Xk : numpy.ndarray
        Matrix with local validation set for segment k.
    PRM : numpy.ndarray
        Projecton matrix for orthogonalization of residuals.

    Returns:
    --------
    numpy.ndarray
        A matrix with orthogonal part for Xpv
    """

    nobj = qk.size

    # project Xk to a set of random vectors and normalize columns
    Z = np.random.randn(nobj, nobj)
    Xpv_orth = np.dot(Z, Xk)
    Xpv_orth = Xpv_orth / np.linalg.norm(Xpv_orth, axis = 0)

    # orthogonalize and rescale rows so sum of their squares is the same as qk
    Xpv_orth = np.dot(Xpv_orth, PRM)
    Xpv_orth = Xpv_orth * (np.sqrt(qk) / np.linalg.norm(Xpv_orth, axis = 1)).reshape(-1, 1)

    return Xpv_orth

def normalize(X: np.ndarray, axis: int = 0) -> np.ndarray:
    """ Normalize rows or columns of X to Euclidean norm. """
    norm = np.sqrt((X * X).sum(axis = axis))
    return X / norm

def simpls(X: np.ndarray, Y: np.ndarray, ncomp: int):
    """
#' SIMPLS algorithm
#'
#' @description
#' SIMPLS algorithm for calibration of PLS model
#'
#' @param X
#' a matrix with x values (predictors)
#' @param Y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#'
#' @return
#' a list with computed weights, x- and y-loadings for PLS regression model.
#'
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression.
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
    """

    (nrows, ncols) = X.shape
    nresp = Y.shape[1]

    S = np.dot(X.T, Y)
    M = np.dot(X.T, X)

    C = np.zeros((nresp, ncomp))
    R = np.zeros((ncols, ncomp))
    P = np.zeros((ncols, ncomp))
    V = np.zeros((ncols, ncomp))
    T = np.zeros((nrows, ncomp))

    for a in range(0, ncomp):
        Vl, sl, Ul = svd(S, full_matrices=False)
        r = Vl[..., :1]
        t = np.dot(X, r)
        tnorm = math.sqrt((t * t).sum())
        t = t / tnorm
        r = r / tnorm

        p = np.dot(X.T, t)
        c = np.dot(Y.T, t)
        u = np.dot(Y, c)
        v = p

        if (a > 0):
            v = v - np.dot(V, np.dot(V.T, p))
            u = u - np.dot(T, np.dot(T.T, u))

        v = v / math.sqrt((v * v).sum())

        R[..., a] = r.flatten()
        V[..., a] = v.flatten()
        P[..., a] = p.flatten()
        T[..., a] = t.flatten()
        C[..., a] = c.flatten()

        M = M - np.dot(p, p.T)
        S = S - np.dot(v, np.dot(v.T, S))

    return(R, P, C)