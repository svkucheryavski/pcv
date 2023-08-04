"""
Helper methods for main procedures.
"""

import numpy as np
import numpy.matlib
import math



def get_splits(nrows: int, nseg: int):
    """
    Generate vector with segment numbers.
    """
    nrep = math.ceil(nrows / nseg)
    return np.matlib.repmat(np.arange(1, nseg + 1, dtype = int), 1, nrep)[0, :nrows]



def get_cvncomp(nrows: int, ncomp: int, nseg: int):
    """
    Return adjusted number of components by taking into account
    size of local calibration sets inside cross-validation loop.
    """
    return min([ncomp, nrows - math.ceil(nrows/nseg) - 1])



def get_cvsettings_manual(cv: np.ndarray, nrows:int, ncomp:int):
    """
    Check the manually entered vector with segment numbers for cross-validation
    and return the list with adjusted CV settings: (indices, nrows, ncomps).
    """

    if cv.size != nrows:
        raise ValueError('if argument "cv" is a vector of indices it must have the same size as value of argument "nrows".')
    nseg = max(cv)
    return (cv, get_cvncomp(nrows, ncomp, nseg), nseg)



def get_cvsettings_rand(nrows:int, ncomp:int, nseg:int):
    """
    Returns a list with main CV settings: (indices, nrows, ncomps) for random splits.
    """

    # shuffle and return the indices and other settings
    cvind = get_splits(nrows, nseg)
    np.random.shuffle(cvind)
    return (cvind, get_cvncomp(nrows, ncomp, nseg), nseg)



def get_cvsettings_ven(nrows:int, ncomp:int, nseg:int, resp:np.ndarray):
    """
    Returns a list with main CV settings: (indices, nrows, ncomps) for Venetian blinds splits.
    """

    if resp is None:
        resp = np.arange(1, nrows + 1)

    if resp.size != nrows:
        raise ValueError('argument "resp" must have the same size as value of argument "nrows".')

    # shuffle and return the indices and other settings
    cvind = get_splits(nrows, nseg)
    row_ind = np.argsort(np.argsort(resp.flatten(), kind = 'mergesort'), kind = 'mergesort')
    return (cvind[row_ind], get_cvncomp(nrows, ncomp, nseg), nseg)



def get_cvsettings(cv, nrows: int, ncomp: int, resp: np.ndarray = None):
    """
    Returns a list with main CV settings: (indices, nrows, ncomps) depending on user input
    and size of calibration set.
    """

    # if user already provided vector with segment numbers - return it
    if isinstance(cv, np.ndarray):
        return get_cvsettings_manual(cv, nrows, ncomp)

    # if not - check that cv is a dictionary and it has at least field "type"
    if type(cv) is not dict:
        raise ValueError('argument "cv" must be either a vector of integers or a dictionary.')

    if 'type' not in cv:
        raise ValueError('if argument "cv" is a dictionary it must have a key "type".')

    # if leave-one-out - return a sequence 1...nrows
    if cv['type'] == 'loo':
        return (np.arange(1, nrows + 1, dtype = int), get_cvncomp(nrows, ncomp, nrows), nrows)

    # if not - first check that the number of segments is provided correctly
    if 'nseg' not in cv or cv['nseg'] < 2 or cv['nseg'] > nrows:
        raise ValueError('wrong value for key "nseg" of argument "cv".')

    # if random splits - return the cv-indices shuffled
    if cv['type'] == 'rand':
        return get_cvsettings_rand(nrows, ncomp, cv['nseg'])

    # if venetian blinds - return the cv-indices ordered based on response values order
    if cv['type'] == 'ven':
        return get_cvsettings_ven(nrows, ncomp, cv['nseg'], resp)


    raise ValueError('wrong value for key "type" of argument "cv".')
