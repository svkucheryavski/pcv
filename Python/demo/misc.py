###################################################################################
#  This file contains several handy functions used in demo examples from demo.py  #
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from prcv.methods import simpls
from scipy.linalg import svd


def scale(X:np.ndarray, center:bool = True, scale:bool = False) -> np.ndarray:
    """ Autoscale columns of X. """

    ncols = X.shape[1]
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    return (X - mX) / sX


def pls_fit(X:np.ndarray, Y:np.ndarray, ncomp:int, center:bool = True, scale:bool = False) -> dict:
    """ Fits PLS model """

    (nrows, ncols) = X.shape

    # correct shape of Y just in case
    if Y.shape is (nrows, ):
        Y.reshape(nrows, 1)

    nresp = Y.shape[1]

    # autoscale data and save vectors for centring and scaling
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    mY = Y.mean(axis = 0) if center else np.zeros(nresp)
    sY = Y.std(axis = 0, ddof = 1) if scale else np.ones(nresp)

    X = (X - mX) / sX
    Y = (Y - mY) / sY

    R, P, C = simpls(X, Y, ncomp)

    return {
        'R': R,
        'P': P,
        'C': C,
        'mX': mX,
        'sX': sX,
        'mY': mY,
        'sY': sY,
        'ncomp': ncomp
    }


def pls_predict(m:dict, X:np.ndarray, Y:np.ndarray) -> dict:
    """ Apply PLS model and compute main outcomes (predictions, errors, variances). """

    nrows, ncols = X.shape
    ncomp = m['ncomp']

    # autoscale predictors
    X = (X - m['mX']) / m['sX']
    T = np.dot(X, m['R'])

    Yp = np.zeros((nrows, ncomp))
    Q = np.zeros((nrows, ncomp))

    # make predictions for a = 1...ncomp components in a model
    for a in range(1, ncomp + 1):
        Pa = m['P'][..., :a]
        Ca = m['C'][..., :a]
        Ta = T[..., :a]
        Ea = X - np.dot(Ta, Pa.T)
        Q[..., a - 1] = (Ea * Ea).sum(axis = 1)
        Yp[..., a - 1] = np.dot(Ta, Ca.T).flatten()


    Yp = Yp * m['sY'] + m['mY']
    SSE = ((Y - Yp) ** 2).sum(axis = 0)
    SSY = ((Y - Y.mean(axis = 0)) ** 2).sum()
    SSX = (X ** 2).sum()

    return {
        'Y': Y,
        'Yp': Yp,
        'RMSE': np.sqrt(SSE / nrows),
        'xexpvar': 1 - Q.sum(axis = 0) / SSX,
        'yexpvar': 1 - SSE / SSY
    }


def pca_fit(X:np.ndarray, ncomp:int, center:bool = True, scale:bool = False) -> dict:
    """ Fits PCA model """

    nrows, ncols = X.shape

    # autoscale data and save vectors for centring and scaling
    mX = X.mean(axis = 0) if center else np.zeros(ncols)
    sX = X.std(axis = 0, ddof = 1) if scale else np.ones(ncols)
    X = (X - mX) / sX

    # get and return PCA outcomes by using SVD
    U, s, V = svd(X, full_matrices=False)
    return {
        'P': V[:ncomp, ...].T,
        'eigenvals': (s[:ncomp] * s[:ncomp]) / (nrows - 1),
        'mX': mX,
        'sX': sX,
        'ncomp': ncomp
    }


def pca_predict(m:dict, X:np.ndarray) -> dict:
    """ Project data to PCA model and return main PCA outcomes (scores, distances, variances). """

    nrows, ncols = X.shape
    ncomp = m['ncomp']

    # autoscale data and save vectors for centring and scaling
    X = (X - m['mX']) / m['sX']

    # get and return PCA outcomes by using SVD
    T = np.dot(X, m['P'])
    U = T / np.sqrt(m['eigenvals'])

    H = np.zeros((nrows, ncomp))
    Q = np.zeros((nrows, ncomp))
    totvar = (X * X).sum()

    for a in range(1, m['ncomp'] + 1):
        Pa = m['P'][..., :a]
        Ta = T[..., :a]
        Ua = U[..., :a]
        Ea = X - np.dot(Ta, Pa.T)
        Q[..., a - 1] = (Ea * Ea).sum(axis = 1)
        H[..., a - 1] = (Ua * Ua).sum(axis = 1)

    return {
        'T': T,
        'H': H,
        'Q': Q,
        'expvar': 1 - Q.sum(axis = 0) / totvar
    }


def plot_distances(res_cal:dict, res_pv:dict, ncomp:int, ax:plt.Axes = None, xlim:[] = None, ylim:[] = None):
    """ Show distance plot (orthogonal vs score) for calibration and PV-set PCA results. """

    if ax is None:
        fig, ax = plt.subplots()

    h_cal = res_cal['H'][..., ncomp - 1]
    q_cal = res_cal['Q'][..., ncomp - 1]

    h_pv = res_pv['H'][..., ncomp - 1]
    q_pv = res_pv['Q'][..., ncomp - 1]

    h0 = h_cal.mean()
    q0 = q_cal.mean()

    ax.grid(True)
    ax.scatter(h_cal / h0, q_cal / q0, label = 'cal', marker = 's', color = 'None', edgecolor = 'blue')
    ax.scatter(h_pv / h0, q_pv / q0, label = 'pv', marker = 'x', color = 'red')

    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    ax.legend()
    ax.set_xlabel('Score distance, h/h0')
    ax.set_ylabel('Orthogonal distance, q/q0')
    ax.set_title('Distance plot (A = ' + str(ncomp) + ')')


def plot_plsres(res_cal:dict, res_pv:dict, ax:plt.Axes = None, type = "xexpvar", title = "X variance"):
    """ Show plot of a PLS result outcome (variance, RMSE, etc.) vs components. """

    if ax is None:
        fig, ax = plt.subplots()

    var_cal = res_cal[type]
    var_pv = res_pv[type]
    comp_seq = range(1, var_cal.size + 1)

    ax.grid(True)
    ax.plot(comp_seq, var_cal, label = 'cal', marker = 's', color = 'blue')
    ax.plot(comp_seq, var_pv, label = 'pv', marker = 'x', color = 'red')

    ax.legend()
    ax.set_xlabel('Components')
    ax.set_ylabel(type)
    ax.set_title(title)


def plot_predictions(res_cal:dict, res_pv:dict, ncomp:int, ax:plt.Axes = None):
    """ Show plot with predicted vs measured values for PLS results. """

    if ax is None:
        fig, ax = plt.subplots()

    y = res_cal['Y']
    yp_cal = res_cal['Yp'][..., ncomp - 1]
    yp_pv = res_pv['Yp'][..., ncomp - 1]

    ax.grid(True)
    ax.scatter(y, yp_cal, label = 'cal', marker = 's', color = 'None', edgecolor = 'blue')
    ax.scatter(y, yp_pv, label = 'pv', marker = 'x', color = 'red')

    ax.legend()
    ax.set_xlabel('y, reference')
    ax.set_ylabel('y, predicted')
    ax.set_title('Predictions (A = ' + str(ncomp) + ')')


