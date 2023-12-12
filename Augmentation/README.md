# Procrustes cross-validation for data augmentation

Our recent pre-print ([arXiv:2312.04911v1](https://arxiv.org/abs/2312.04911)) demonstrates that PCV, in addition to its original purpose, can also address the problem of data augmentation. The original paper is now submitted to one of the peer-review journals, we will update this text when it is accepted.

Meanwhile this part of repository contains two Jupyter notebooks that reproduce some of the examples shown in the pre-print in a simplified way.

## How to use

File `Tector.ipynb` contains notebook with code reproducing the examples with *Tecator* dataset. It trains 2 x 10 simple ANN regression models and applies each model to independent test set. Half of the models are trained using the original calibration set (170 objects). Second half — the calibration set augmented with 10 PV-sets generated using Procrustes cross-validation with radom splits to 10 segments and PLS-decomposition with 10 latent variables. The code saves RMSEP values obtained for each model and shows box plot for the values at the end. To run the code you need to download `Tecator.csv` file with dataset (also located here).

Second notebook, `Heart.ipynb`, reproduces some of the examples with *Heart* dataset in the same way. To use this notebook you need to download `Heart.csv` file with dataset, where all categorical variables already converted to dummy. The notebook implements PLS-version of the augmentation algorithm (in the pre-print we use both PLS and PCA versions for this example) and uses the same train/test splitting strategy as described in the pre-print.

Code works best if you have [CUDA device](https://en.wikipedia.org/wiki/CUDA), it will work without the device as well, but calculations will take longer time (on MacBook Air M2 it takes about 10-12 minutes to run all examples from the Tecator notebook).

The Tecator dataset was originally taken by Tecator Infratec Food and Feed Analyzer working by the Near Infrared Transmission (NIT) principle. The dataset was downloaded from the StatLib public dataset archive (http://lib.stat.cmu.edu/datasets/). The Heart dataset was downloaded from UC Irvine Machine Learning Repository and adapted for the purpose of the paper. More details can be found in the pre-print.

## Short description

The general idea is simple. If we use random segmented cross-validation, we can obtain a large number of PV-sets for the same conditions. They will all be different but will capture the variance-covariance structure (full or partial) of the original training set, similar to having many samples taken from the same population. What if we simply combine them all with a training set, thereby augmenting the original training data?

Results show that, for example, in the case of Tecator data and a straightforward ANN model with several layers, using data augmented this way significantly reduces the RMSE computed for the independent test set (in some scenarios, up to a 2-3 times reduction). This holds true for most of the spectroscopic examples I have explored so far. Sometimes the reduction is relatively small (10-20%), sometimes quite large (several times). Experiments with tabular data also demonstrate improvement, such as in discrimination accuracy.

We can also observe that for models robust to overfitting, such as Random Forest, augmentation does not yield any improvements. Overall, there is still room for investigation and potential improvements. However, what we have seen so far fills us with optimism, and we encourage everyone to try how the method works on your data, especially if you have a dataset with a moderate to high degree of collinearity, a relatively small number of data points (from tens to a few hundreds), and you want to use ANN for modeling.


