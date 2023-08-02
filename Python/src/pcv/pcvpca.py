import numpy as np
import numpy.linalg
import math
from .misc import get_xpvorth, get_cvsettings

def get_pcamodel(X: np.ndarray, ncomp: int):
    """
    Create PCA model by using SVD and truncating loadings and singular values.
    """
    U, s, V = np.linalg.svd(X)
    P = V.transpose()[..., :ncomp]
    s = s[:ncomp]

    return (P, s)



def pcvpca(X: np.ndarray,
        ncomp: int,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        cvscope: str = 'global'):
    """
    Procrustes cross-validation based on PCA decomposition.

    Parameters:
    -----------
    X : numpy.ndarray
        Global calibration set as 2D Numpy array.
    ncomp : int
        Number of components to use for generation (must be larger than expected optimal number).
    center : bool
        Logical, center or not columns of X (better do it).
    scale : bool
        Logical, standardize or not columns of X.
    cv : array|dict
        Parameters of cross-validation loop, see details below.
    cvscope: str
        Defines how local calibration set must be centered and/or scaled.

    Returns:
    --------
    Xpv : numpy.ndarray
        2D array of the same size as X with generated PV-set.


    """

    (nrows, ncols) = X.shape

    # compute global mean and standard deviation and autoscale the whole data
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    X = (X - mX) / sX

    # indices for cross-validation
    (cvind, cvncomp, cvnseg) = get_cvsettings(cv, nrows, ncomp)

    # create a global model
    P, s = get_pcamodel(X, cvncomp)
    PRM = np.eye(ncols) - np.dot(P, np.transpose(P))

    # prepare empty matrix for pseudo-validation set
    Xpv = np.zeros((nrows, ncols))

    for k in range(1, cvnseg + 1):

        indc = cvind != k
        indk = cvind == k

        Xc = X[indc, ...]
        Xk = X[indk, ...]

        if cvscope == 'local':
            mXl = Xc.mean(axis = 0) if center else np.zeros(ncols)
            sXl = Xc.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
            Xc = (Xc - mXl) / sXl
            Xk = (Xk - mXl) / sXl

        # get local model and correct direction of the loadings
        Pk, sk = get_pcamodel(Xc, cvncomp)
        a = np.arccos((P * Pk).sum(axis = 0)) < math.pi / 2
        Pk = Pk * (2 * a - 1)

        # get scores and residuals by projection local validation set to the local model
        Tk = np.dot(Xk, Pk)
        Ek = Xk - np.dot(Tk, np.transpose(Pk))
        qk = (Ek * Ek).sum(axis = 1)

        # compute the parallel component of Xpv
        Xpv_hat = np.dot(Tk, np.transpose(P))

        # compute the orthogonal component of Xpv
        Xpv_orth = get_xpvorth(qk, Xk, PRM)

        # create and save the Xpv
        Xpv[indk, ...] = Xpv_hat + Xpv_orth

    # uscenter and unscale the data using global mean and std and return the result
    return Xpv * sX + mX


