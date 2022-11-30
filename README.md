# Procrustes cross-validation


Procrustes cross-validation is a new approach for validation of chemometric models. It makes possible to generate a new dataset, named *PV-set* and use it for validation of models in the same way as with an independent validation set.

This repository contains implementation of Procrustes cross-validation for R and MATLAB. There are two versions of the method:

1. The *original* version was created in 2020. This repository contains both code and references for the original version, so if you used it before you can safely continue. No breaking changes are introduced to this code.

2. The *new* version was created in 2022. The code for this version is already available in this repository, however the paper, which describes the idea behind the new version, is currently under peer-review. We will update this text as soon as the paper is out.

All details for both versions are given below.

## The new version

The new version of the method is a natural development of the original idea. However, we generalized it, which made the implementation much faster for PCA/SIMCA models (several order of magnitude) and allowed us to extend the method for validation of regression models (PCR, PLS1 and PLS2).

Moreover, the implementation of the new version (R and MATLAB code you will find here) is more flexible in terms of splitting the data into segments, compared to the implementation we made for the original version. Now you can chose between leave-one-out, systematic splits, random splits as well as provide your own vector with segment indices.

This version is implemented as a set of the following methods (both in R and MATLAB): `pcvpca()`, `pcvpcr()` and `pcvpls()`. The first implementation is for PCA/SIMCA models, the second is for Principal Component Regression and the last one is for Partial Least Squares regression (both for single and multivariate response).

### Reference

The paper for the new version of the method is currently under peer-review. We will update this information as soon as it is out:

Kucheryavskiy, S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation of multivariate regression model*. Submitted, 2022.

### How to use in R

All code is located in `pcv.R`, which you will find inside folder `R` in this repository, no extra packages are needed. You can also download a zip archive with both R and MATLAB versions from the Release section. The simply copy the `pcv.R` file to your working directory and use `source("pcv.R")` inside your script to load all necessary functions to the environment.

Below are examples of syntax with all available parameters:

```r
# for PCA/SIMCA models
Xpv <- pcvpca(X, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PCR models
Xpv <- pcvpcr(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PLS models
Xpv <- pcvpls(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix, not a data frame). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. The method generates PV-set only for predictors, the response values for the calibration set and the PV-set are the same.

Parameter `ncomp` is a number of principal components in case of PCA/PCR models or number of latent variables in a PLS model. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `center` and `scale` define if the predictors must be mean centered and/or standardized. By default `center = TRUE` and `scale = FALSE`. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with the calibration set.

Finally, parameter `cv` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `cv = 4`). In this case this number specifies number of segments for random splits, except `cv = 1` which is a special case for leave-one-out (full cross-validation).

* A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make the split, you can select one of the following: `"full"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in the code examples above, tells PCV to use Venetian blinds splits with 4 segments.

* A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number of values in this vector must be the same as number of rows in the training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.


File `demo.R` contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.

## How to use in MATLAB

All code is located in `MATLAB` folder, no extra packages are needed. Download this folder to your computer (you can get it as ZIP file in Release section), unzip and rename as you want. After that run MATLAB, find the folder in file browser inside MATLAB, click right mouse button on that folder and select "Add ...".

Below you will see examples of syntax with all parameters available:

```matlab
% for PCA/SIMCA models
Xpv = pcvpca(X, 20, 'Center', true, 'Scale', false, 'CV', {"ven", 4});

% for PCR models
Xpv = pcvpcr(X, y, 20, 'Center', true, 'Scale', false, 'CV', {"ven", 4});

% for PLS models
Xpv = pcvpls(X, y, 20, 'Center', true, 'Scale', false, 'CV', {"ven", 4});
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. The method generates PV-set only for predictors, the response values for calibration set and PV-set are the same.

The second parameter is a number of principal components in case of PCA/PCR models or number of latent variables in PLS model. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `Center` and `Scale` defines if the predictors must be mean centered and/or standardized. By default `Center` is set to `true` and `Scale` is set to  `false`. Regardless which settings you use the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with calibration set.

Finally, parameter `CV` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `'CV', 4`). In this case this number specifies number of segments for random splits, except `'CV', 1` which is a special case for leave-one-out (full cross-validation).

* A cell vector with 2 values: `{"name", nseg}`. In this case `"name"` defines the way to make the split, you can select one of the following: `"full"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into.

* A vector with integer numbers, e.g. `'CV', [1, 2, 3, 1, 2, 3, 1, 2, 3]`. In this case number of values in this vector must be the same as number of rows in training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.


File `demo.m` contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.


## The original (old) version

The original version is available is a function `pcv`. It will remain unchanged to ensure compatibility, so if you have already used the original version, you can safely continue. However if you start a new project we encourage you to switch to the new version of the algorithm.

### Reference

The original version has been published in this paper:

Kucheryavskiy, S., Zhilin, S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation—A Bridge between Cross-Validation and Independent Validation Sets.* Analytical Chemistry,  92 (17), 2020. pp.11842–11850. DOI: [10.1021/acs.analchem.0c02175](https://doi.org/10.1021/acs.analchem.0c02175)

Please cite this paper only if indeed you use the original version of the algorithm, otherwise, use the reference for the new version (see above).

The next paper shows applicability and efficiency of the proposed method in case of short datasets:

Pomerantsev A. L., Rodionova O. Ye. *Procrustes Cross-Validation of short datasets in PCA context.*
Talanta, 226, 2021. DOI: [10.1016/j.talanta.2021.122104](https://doi.org/10.1016/j.talanta.2021.122104)

### How to use in R

All code is located in `pcvold.R`, no extra packages are needed. Simply copy the file to your working directory and use `source("pcvold.R")` inside your script to load all necessary functions to the environment. Then the syntax is follows:

```r
pvset <- pcv(X, ncomp, nseg, scale)
```

Here `X` is your calibration set (as a numerical matrix), `ncomp` is a number of principal components for PCA decomposition (use at lease enough to explain 95-99% of variation of the data values), `nseg` is the number of segments for cross-validation procedure. So far only systematic cross-validation (venetian blinds) is supported, so make sure that rows of `X` are sorted correctly or shuffled. Parameter `scale` allows to standardize your data prior to the generation, which is useful if your variables have different nature and/or units. The generated data will be unscaled though.

File `demoold.R` contains a demo code based on *NIRSim* dataset from the original paper. See comments in the code for more details.

### How to use in MATLAB

See details how to get the MATLAB version of the new version. The folder contains implementation of the original version as well. Use the following syntax:

```matlab
pvset = pcv(X, nComp, nSeg, Scale);
```

File `demoold.m` contains a demo code based on *NIRSim* dataset from the paper. See comments in the code for more details.


## Bugs and improvements

The code will be improved and extended gradually. If you found a bug please report using issues or send an email.
