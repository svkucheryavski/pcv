#####################################################################################
#                                                                                   #
#  This file contains demo code for implementation of Procrustes cross-validation   #
#  method in Python. The code creates plots similar to what you can see on figures  #
#  with Corn example in this paper: https://doi.org/10.1016/j.aca.2023.341096       #
#                                                                                   #
#  Check this for more details: https://github.com/svkucheryavski/pcv/        #
#                                                                                   #
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt

from prcv.methods import pcvpca, pcvpcr, pcvpls

# load helper functions from misc.py (you need to download it as well)
from misc import scale, pca_fit, pca_predict, plot_distances, \
    pls_fit, pls_predict, plot_plsres, plot_predictions


# load original spectra and response values for Corn dataset
# remember to download the CSV file from the repo
D = np.genfromtxt('corn.csv', delimiter=',')
X = D[:, 1:]
Y = D[:, :1]

# 1. PCA based examples

## create pseudo-validation set using PCA based algorithm
Xpv = pcvpca(X, ncomp = 30, cv = {'type': 'ven', 'nseg': 4})

## make PCA model for calibration set
m = pca_fit(X, 20)

## project calibration and PV-set to the PCA model
res_cal = pca_predict(m, X)
res_pv = pca_predict(m, Xpv)

## show distance plots for A = 2 and A = 20
## similar to ones shown in Supplementary materials
fig, ((ax1, ax2)) = plt.subplots(nrows = 1, ncols = 2)
plot_distances(res_cal = res_cal, res_pv = res_pv, ncomp =  2, xlim = [0, 22], ylim = [0, 12], ax = ax1)
plot_distances(res_cal = res_cal, res_pv = res_pv, ncomp = 20, xlim = [0,  7], ylim = [0, 12], ax = ax2)
plt.show()


# 2. PCR based examples

## create pseudo-validation set using PCR based algorithm
Xpv, D = pcvpcr(X, Y, ncomp = 20, cv = {'type': 'ven', 'nseg': 4})

## show plot with original and generated spectra
## similar to plots shown in Figure 2

plt.subplot(121)
plt.plot(scale(X).T)
plt.title('Calibration set (mean centered)')

plt.subplot(122)
plt.plot(scale(Xpv).T)
plt.title('PV-set (mean centered)')

plt.show()


# 3. PLS based examples

## create pseudo-validation set for PLS model
Xpv, D = pcvpls(X, Y, ncomp = 20, cv = {'type': 'ven', 'nseg': 4})

## create PLS model
m = pls_fit(X, Y, ncomp = 20)

## apply PLS model to calibration and PV-set
res_cal = pls_predict(m, X, Y)
res_pv = pls_predict(m, Xpv, Y)

## show RMSE and predicted vs. measured plots
## similar to plots shown in Figure 5
fig, ((ax1, ax2, ax3)) = plt.subplots(nrows = 1, ncols = 3)
plot_plsres(res_cal = res_cal, res_pv = res_pv, ax = ax1, type = "RMSE", title = "RMSE")
plot_predictions(res_cal = res_cal, res_pv = res_pv, ax = ax2, ncomp = 10)
plot_predictions(res_cal = res_cal, res_pv = res_pv, ax = ax3, ncomp = 20)
plt.show()


## show heatmap and boxplot for elements of D
## similar to the first plot shown in Figure 6

plt.subplot(211)
plt.imshow(D[::-1, :], aspect = 'auto', extent = (0.5, D.shape[1] + 0.5, 0.5, D.shape[0] + 0.5))
plt.xticks([])
plt.yticks(range(1, D.shape[0] + 1))
plt.ylabel('Segments, k')

plt.subplot(212)
plt.boxplot(D)
plt.ylabel('c_k/c')
plt.xlabel('Components, a')

plt.show()
