"""
Procrustes cross-validation for three-way data: FU, BU and HSV methods.
"""

import contextlib
import warnings
from typing import Optional

import numpy as np
from scipy.linalg import svd
from scipy.optimize import linear_sum_assignment

from .methods import pcvpca
from .misc import get_cvsettings


_EPS = np.finfo(float).eps



################################
# Main methods                 #
################################

def pcv3wfu(X: np.ndarray,
        ncomp: Optional[int] = None,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        seed: Optional[int] = None) -> np.ndarray:
    """
    Procrustes cross-validation for three-way data based on full unfolding (FU).

    Parameters:
    -----------
    X : numpy.ndarray
        Calibration set as 3D array with shape (I, J, K): samples x mode 2 x mode 3.
    ncomp : int, default = min(I - 1, J * K, 30)
        Number of PCA components (more than the expected optimal number of components of the
        three-way model).
    center : bool, default = True
        Logical, center the unfolded data or not.
    scale : bool, default = False
        Logical, scale the unfolded data or not (each of the J * K columns separately).
    cv : array|dict, default = {'type': 'ven', 'nseg': 4}
        Parameters of cross-validation loop, see `pcvpca()` for details.
    seed : int, default = None
        Makes random splits and random orthogonal parts reproducible. The random state of
        NumPy is restored afterwards.

    Returns:
    --------
    Xpv : numpy.ndarray
        3D array of the same shape as X with generated PV-set.

    Notes:
    ------
    The method unfolds X to a matrix with I rows and J * K columns, generates PV-set for the
    matrix using `pcvpca()` with global CV scope, and folds the result back to a 3D array.
    The PV-set reproduces cross-validation of the PCA model of the unfolded data. Use enough
    components, so the PCA subspace contains the subspace of the three-way models (e.g. PARAFAC
    or Tucker3) you are going to validate.
    """

    (nrows, n2, n3) = _getdim3w(X)
    if ncomp is None:
        ncomp = min(nrows - 1, n2 * n3, 30)

    with _seed3w(seed):
        cvind = _splitcv3w(cv, nrows)
        Xpv = pcvpca(X.reshape(nrows, n2 * n3), ncomp, center = center, scale = scale, cv = cvind,
            cvscope = 'global')

    return Xpv.reshape(nrows, n2, n3)



def pcv3wbu(X: np.ndarray,
        mode: int = 2,
        ncomp: Optional[int] = None,
        center: bool = True,
        scale: bool = False,
        cv: dict = {'type': 'ven', 'nseg': 4},
        seed: Optional[int] = None) -> np.ndarray:
    """
    Procrustes cross-validation for three-way data based on block unfolding (BU).

    Parameters:
    -----------
    X : numpy.ndarray
        Calibration set as 3D array with shape (I, J, K): samples x mode 2 x mode 3.
    mode : int, default = 2
        Which mode to split into blocks (2 or 3).
    ncomp : int, default = min(I - 1, size of block, 30)
        Number of PCA components for each block.
    center : bool, default = True
        Logical, center each block or not.
    scale : bool, default = False
        Logical, scale each block or not.
    cv : array|dict, default = {'type': 'ven', 'nseg': 4}
        Parameters of cross-validation loop, see `pcvpca()` for details.
    seed : int, default = None
        Makes random splits and random orthogonal parts reproducible. The random state of
        NumPy is restored afterwards.

    Returns:
    --------
    Xpv : numpy.ndarray
        3D array of the same shape as X with generated PV-set.

    Notes:
    ------
    The method splits X into J matrices of size I x K (if `mode = 2`) or K matrices of size
    I x J (if `mode = 3`). Then it generates PV-set for each matrix using `pcvpca()` with global
    CV scope and combines the PV-sets into a 3D array. All blocks use the same split of the
    samples into segments.
    """

    (nrows, n2, n3) = _getdim3w(X)
    if mode not in (2, 3):
        raise ValueError('argument "mode" must be 2 or 3.')

    if ncomp is None:
        ncomp = min(nrows - 1, n3 if mode == 2 else n2, 30)

    Xpv = np.zeros((nrows, n2, n3))
    with _seed3w(seed):
        cvind = _splitcv3w(cv, nrows)
        for b in range(n2 if mode == 2 else n3):
            Xb = X[:, b, :] if mode == 2 else X[:, :, b]
            Xpvb = pcvpca(Xb, ncomp, center = center, scale = scale, cv = cvind, cvscope = 'global')
            if mode == 2:
                Xpv[:, b, :] = Xpvb
            else:
                Xpv[:, :, b] = Xpvb

    return Xpv



def pcv3whsv(X: np.ndarray,
        ncomp: tuple,
        center: bool = True,
        scale = 'none',
        cv: dict = {'type': 'ven', 'nseg': 4},
        orth: str = 'random',
        align: str = 'optimal',
        seed: Optional[int] = None):
    """
    Procrustes cross-validation for three-way data based on HO-SVD (HSV).

    Parameters:
    -----------
    X : numpy.ndarray
        Calibration set as 3D array with shape (I, J, K): samples x mode 2 x mode 3. All values
        must be finite.
    ncomp : tuple
        Two positive integers, (R2, R3), number of components in modes 2 and 3.
    center : bool, default = True
        Logical, subtract the mean slab or not.
    scale : bool|str, default = 'none'
        Scaling within modes: 'none' (or False), 'mode2', 'mode3' or 'both' (or True). A variable
        with zero variance (also when it is constant over the samples) cannot be scaled and gives
        an error.
    cv : array|dict, default = {'type': 'ven', 'nseg': 4}
        Parameters of cross-validation loop, see `pcvpca()` for details.
    orth : str, default = 'random'
        How to construct the orthogonal part: 'random' (random mix of the segment's slabs, similar
        to the two-way method) or 'canonical' (each slab's own local residual, no randomness).
    align : str, default = 'optimal'
        How to match local components to the global ones: 'optimal' (maximal sum of absolute
        cosines) or 'greedy' (global components in order of variance choose first).
    seed : int, default = None
        Makes random splits and random orthogonal parts reproducible. The random state of
        NumPy is restored afterwards.

    Returns:
    --------
    Xpv : numpy.ndarray
        3D array of the same shape as X with generated PV-set.
    info : dict
        Details: 'ncomp' (numbers of components actually used), 'cvind' (segment number of each
        sample), 'orth', 'align', 'scale', 'prep' (global preprocessing parameters), 'alignment'
        and 'valid_subranks' (see notes).

    Notes:
    ------
    The method models X with truncated HO-SVD in modes 2 and 3 (a Tucker2 model) and generates
    PV-set, which, when projected on the global model, gives the same scores and orthogonal
    distances as cross-validation of the model, sample by sample.

    Scaling divides mode 2 and/or mode 3 variables by their standard deviation after centering.
    With 'both', mode 2 is scaled first and mode 3 is then scaled using the mode 2 scaled data.
    The PV-set is returned in the original units.

    `info['alignment']` is a 2D array with columns: segment, mode, global component, matched
    local component and absolute cosine between the two. Segments and components are numbered
    from 1.

    A PV-set generated with (R2, R3) components also reproduces cross-validation of the models
    with fewer components, r2 <= R2 and r3 <= R3, unless the local components are matched to the
    global ones across the boundary of the first r2 (or r3) components. The numbers of components
    for which this holds are given in `info['valid_subranks']`, and the method gives a warning if
    some of them are missing.
    """

    if orth not in ('random', 'canonical'):
        raise ValueError('argument "orth" must be "random" or "canonical".')

    if align not in ('optimal', 'greedy'):
        raise ValueError('argument "align" must be "optimal" or "greedy".')

    scale = _normalizescale3w(scale)
    _getdim3w(X)

    if not np.all(np.isfinite(X)):
        raise ValueError('argument "X" contains NaN or infinite values, interpolate or impute them first.')

    ncomp = _checkncomp3w(ncomp)

    with _seed3w(seed):
        cvind = _splitcv3w(cv, X.shape[0])
        return _pcv3whsv_run(X.astype(float), ncomp, cvind, center, scale, orth, align)



def _pcv3whsv_run(X: np.ndarray, ncomp: tuple, cvind: np.ndarray, center: bool, scale: str,
        orth: str, align: str):
    """
    Main loop of pcv3whsv() (arguments already checked).
    """

    nseg = int(cvind.max())

    # global preprocessing and global factors
    prep = _prep3w(X, center, scale)
    Xp = prep['Xp']
    (U2g, U3g, r2min, r3min) = _modefactors3w(Xp, ncomp[0], ncomp[1])

    # local factors for every segment
    local = []
    for k in range(1, nseg + 1):
        (U2l, U3l, r2, r3) = _modefactors3w(Xp[cvind != k], ncomp[0], ncomp[1])
        local.append((U2l, U3l))
        r2min = min(r2min, r2)
        r3min = min(r3min, r3)

    # ranks attainable by the global and all local models
    R2 = min(ncomp[0], r2min)
    R3 = min(ncomp[1], r3min)

    if R2 < 1 or R3 < 1:
        raise ValueError('the data cannot supply at least one component in every mode for all segments.')

    if R2 < ncomp[0] or R3 < ncomp[1]:
        warnings.warn(f'pcv3whsv: "ncomp" reduced to ({R2}, {R3}) (maximum for all segments).')

    U2 = U2g[:, :R2]
    U3 = U3g[:, :R3]
    cdir = _complementdir3w(U2, U3)

    Xpv = np.zeros(X.shape)
    alignment = []

    for k in range(1, nseg + 1):
        indk = np.flatnonzero(cvind == k)
        Xk = [Xp[i] for i in indk]

        # local factors matched to the global ones
        (U2k, perm2, cos2) = _align3w(U2, local[k - 1][0][:, :R2], align)
        (U3k, perm3, cos3) = _align3w(U3, local[k - 1][1][:, :R3], align)
        alignment.append(_aligndiag3w(k, 2, perm2, cos2))
        alignment.append(_aligndiag3w(k, 3, perm3, cos3))

        # explained part: local scores with global factors; residuals of the local model
        Xhat = []
        Eres = []
        for Xi in Xk:
            Ti = U2k.T @ Xi @ U3k
            Eres.append(Xi - U2k @ Ti @ U3k.T)
            Xhat.append(U2 @ Ti @ U3.T)

        Eo = _getxpvorth3w(Xk, Eres, U2, U3, cdir, orth)
        for m, i in enumerate(indk):
            Xpv[i] = Xhat[m] + Eo[m]

    Xpv = _prep3w_remove(Xpv, prep)

    # sub-rank r is valid if in every segment components 1..r are matched within 1..r
    alignment = np.vstack(alignment)
    valid2 = _validsubranks3w(alignment, 2, R2)
    valid3 = _validsubranks3w(alignment, 3, R3)

    if len(valid2) < R2 or len(valid3) < R3:
        warnings.warn('pcv3whsv: PV-set does not reproduce cross-validation for some smaller numbers '
            'of components, see info["valid_subranks"].')

    info = {
        'ncomp': (R2, R3),
        'cvind': cvind,
        'orth': orth,
        'align': align,
        'scale': scale,
        'prep': {key: prep[key] for key in ('M', 'S2', 'S3', 'center', 'scale')},
        'alignment': alignment,
        'valid_subranks': {'r2': valid2, 'r3': valid3}
    }

    return (Xpv, info)



################################
# Checks and settings          #
################################

def _getdim3w(X):
    """
    Check that X is a numeric 3D array and return its shape.
    """

    if not isinstance(X, np.ndarray) or X.ndim != 3 or not np.issubdtype(X.dtype, np.number):
        raise ValueError('argument "X" must be a numeric 3D array.')

    if min(X.shape) < 2:
        raise ValueError('argument "X" must have at least 2 values in each mode.')

    return X.shape



def _normalizescale3w(scale) -> str:
    """
    Convert argument "scale" to one of 'none', 'mode2', 'mode3', 'both'.
    """

    if isinstance(scale, (bool, np.bool_)):
        return 'both' if scale else 'none'

    if scale not in ('none', 'mode2', 'mode3', 'both'):
        raise ValueError('argument "scale" must be False, True or one of "none", "mode2", "mode3", "both".')

    return scale



def _checkncomp3w(ncomp) -> tuple:
    """
    Check that ncomp contains two positive integers and return them as a tuple of ints.
    """

    try:
        values = np.asarray(ncomp, dtype = float).ravel()
    except (TypeError, ValueError):
        values = np.array([])

    if values.size != 2 or not np.all(np.isfinite(values)) or np.any(values < 1) or \
            np.any(values != np.round(values)):
        raise ValueError('argument "ncomp" must be two positive integers, (R2, R3).')

    return (int(values[0]), int(values[1]))



def _splitcv3w(cv, nrows: int) -> np.ndarray:
    """
    Segment number of each sample, checked to be 1, 2, ..., nseg.
    """

    if isinstance(cv, np.ndarray):
        cv = cv.ravel()

    (cvind, _, _) = get_cvsettings(cv, nrows, 1)
    cvind = np.asarray(cvind, dtype = float)
    nseg = cvind.max() if cvind.size > 0 and np.all(np.isfinite(cvind)) else np.nan

    if not np.isfinite(nseg) or np.any(cvind != np.round(cvind)) or nseg < 2 or np.any(cvind < 1) or \
            np.unique(cvind).size != nseg:
        raise ValueError('segment numbers in "cv" must be 1, 2, ..., n (n >= 2), each used at least once.')

    return cvind.astype(int)



@contextlib.contextmanager
def _seed3w(seed):
    """
    Set NumPy random seed if it is provided and restore the previous random state afterwards.
    """

    if seed is None:
        yield
        return

    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)



################################
# Preprocessing                #
################################

def _prep3w(X: np.ndarray, center: bool, scale: str) -> dict:
    """
    Center (mean slab) and scale (within modes) a 3D array.
    """

    (_, n2, n3) = X.shape
    M = X.mean(axis = 0) if center else np.zeros((n2, n3))
    S2 = np.ones(n2)
    S3 = np.ones(n3)
    Xp = X - M

    # zero variance check: sd below the rounding noise of centering (10 * I * eps times the
    # largest absolute value of the variable before centering) counts as zero
    A = np.abs(X) * (10 * X.shape[0] * _EPS)

    if scale in ('mode2', 'both'):
        S2 = Xp.transpose(1, 0, 2).reshape(n2, -1).std(axis = 1, ddof = 1)
        _checkstd3w(S2, A.max(axis = (0, 2)), 2)
        Xp = Xp / S2[None, :, None]
        A = A / S2[None, :, None]

    if scale in ('mode3', 'both'):
        S3 = Xp.transpose(2, 0, 1).reshape(n3, -1).std(axis = 1, ddof = 1)
        _checkstd3w(S3, A.max(axis = (0, 1)), 3)
        Xp = Xp / S3[None, None, :]

    return {'Xp': Xp, 'M': M, 'S2': S2, 'S3': S3, 'center': center, 'scale': scale}



def _checkstd3w(s: np.ndarray, tol: np.ndarray, mode: int):
    """
    Raise an error if some variables of a mode have zero variance (sd not above tol).
    """

    bad = np.flatnonzero(~(np.isfinite(s) & (s > tol)))
    if bad.size > 0:
        raise ValueError(f'cannot scale mode {mode}, zero variance in variable(s) ' +
            ', '.join(str(b + 1) for b in bad[:5]) + '.')



def _prep3w_remove(Xp: np.ndarray, prep: dict) -> np.ndarray:
    """
    Undo preprocessing made by _prep3w() (back to original units).
    """
    return Xp * prep['S2'][None, :, None] * prep['S3'][None, None, :] + prep['M']



################################
# Model helpers                #
################################

def _numrank3w(s: np.ndarray, tol: float = np.sqrt(_EPS)) -> int:
    """
    Number of singular values above the relative tolerance.
    """

    if s.size == 0 or s[0] <= 0:
        return 0

    return int(np.sum(s > tol * s[0]))



def _modefactors3w(X: np.ndarray, R2: int, R3: int):
    """
    Mode 2 and mode 3 factors (truncated SVD of unfoldings) and numerical ranks.
    """

    (_, n2, n3) = X.shape
    U2, s2, _ = svd(X.transpose(1, 0, 2).reshape(n2, -1), full_matrices = False)
    U3, s3, _ = svd(X.transpose(2, 0, 1).reshape(n3, -1), full_matrices = False)
    return (U2[:, :R2], U3[:, :R3], _numrank3w(s2), _numrank3w(s3))



def _proj3w(S: np.ndarray, U2: np.ndarray, U3: np.ndarray) -> np.ndarray:
    """
    Projection of a slab onto range(U2) x range(U3).
    """
    return U2 @ (U2.T @ S @ U3) @ U3.T



################################
# Orthogonal part              #
################################

def _complvec3w(U: np.ndarray) -> np.ndarray:
    """
    Unit vector orthogonal to the columns of U.
    """

    j = np.argmin((U * U).sum(axis = 1))
    e = np.zeros((U.shape[0], 1))
    e[j] = 1.0
    w = e - U @ (U.T @ e)
    return w / np.sqrt((w * w).sum())



def _complementdir3w(U2: np.ndarray, U3: np.ndarray):
    """
    Fixed unit direction outside the model space (None if there is none).
    """

    # first factor with sign fixed, so the direction does not depend on the SVD sign convention
    first = lambda U: U[:, :1] * np.sign(U[np.argmax(np.abs(U[:, 0])), 0])

    if U2.shape[1] < U2.shape[0]:
        return _complvec3w(U2) @ first(U3).T

    if U3.shape[1] < U3.shape[0]:
        return first(U2) @ _complvec3w(U3).T

    return None



def _getxpvorth3w(Xk: list, Eres: list, U2: np.ndarray, U3: np.ndarray, cdir, orth: str) -> list:
    """
    Orthogonal parts of the PV-set slabs of one segment.

    Each orthogonal part lies outside the global model space and has the same norm as the local
    residual of the slab. Its direction comes from a random mix of the segment's slabs ('random')
    or from the slab's own local residual ('canonical').
    """

    nk = len(Xk)
    shape = (U2.shape[0], U3.shape[0])
    eps = 100 * _EPS

    Xs = np.stack(Xk)
    Z = np.random.uniform(-1.0, 1.0, (nk, nk)) if orth == 'random' else np.zeros((nk, nk))
    E = []

    for i in range(nk):
        qi = (Eres[i] * Eres[i]).sum()
        if qi <= 0:
            E.append(np.zeros(shape))
            continue

        if orth == 'canonical':
            src = Eres[i]
        else:
            src = np.tensordot(Z[i], Xs, axes = 1)

        Eo = src - _proj3w(src, U2, U3)

        # source lies (almost) inside the model space: use the fixed direction
        if (Eo * Eo).sum() <= eps * (src * src).sum():
            if cdir is None:
                E.append(np.zeros(shape))
                continue
            Eo = cdir

        E.append(Eo * np.sqrt(qi / (Eo * Eo).sum()))

    return E



################################
# Matching components          #
################################

def _align3w(Ug: np.ndarray, Ul: np.ndarray, method: str):
    """
    Match local factors to global ones by permutation and sign.

    Returns the matched local factors, the matched local component for each global one
    (numbered from 1) and the absolute cosines of the matched pairs.
    """

    R = Ug.shape[1]
    if Ul.shape[1] != R:
        raise ValueError('factor matrices must have the same number of columns.')

    M = np.abs(Ug.T @ Ul)
    if method == 'optimal':
        (_, assigned) = linear_sum_assignment(M, maximize = True)
    else:
        assigned = np.zeros(R, dtype = int)
        used = np.zeros(R, dtype = bool)
        for g in range(R):
            cand = np.flatnonzero(~used)
            best = cand[np.argmax(M[g, cand])]
            assigned[g] = best
            used[best] = True

    Ual = Ul[:, assigned]
    sgn = np.sign((Ug * Ual).sum(axis = 0))
    sgn[sgn == 0] = 1

    return (Ual * sgn, assigned + 1, M[np.arange(R), assigned])



def _aligndiag3w(segment: int, mode: int, perm: np.ndarray, cos: np.ndarray) -> np.ndarray:
    """
    Matching diagnostics of one segment and mode: segment, mode, component, local component, cosine.
    """
    R = perm.size
    return np.column_stack((np.full(R, segment), np.full(R, mode), np.arange(1, R + 1), perm, cos))



def _validsubranks3w(alignment: np.ndarray, mode: int, R: int) -> list:
    """
    Numbers of components r for which components 1..r are matched within 1..r in every segment.
    """

    ok = np.ones(R, dtype = bool)
    rows = alignment[alignment[:, 1] == mode]
    for k in np.unique(rows[:, 0]):
        perm = rows[rows[:, 0] == k, 3]
        ok &= np.maximum.accumulate(perm) == np.arange(1, R + 1)

    return [int(r) for r in np.flatnonzero(ok) + 1]
