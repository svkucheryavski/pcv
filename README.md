# Procrustes cross-validation

<img src="logo.png" width="75" height="75" style="float:left;padding-right:10px;">Procrustes cross-validation (PCV) is a new approach for validation of chemometric models. It generates a new dataset, named *PV-set*, which can be used for validation of PCA/SIMCA/PCR/PLS models in the same way as with an independent validation set.

This repository contains source code with implementation of Procrustes cross-validation for R and MATLAB as well as some practical details. There are two versions of the method:

1. The *original* version was created in 2020. It is based on rotations in variable space, is a bit slow, and works only for PCA/SIMCA models. Paper [1] in References section below describes the original method and paper [2] shows some practical examples of its use for SIMCA classification.

2. The *new* version was created in 2022 as a generalization of the original approach. The new version is much faster and can also work with regression models (e.g. PCR, PLS1, PLS2). Paper [3] describes the new method, however currently the paper is under peer-review in a journal, we will update the information when it is out.

Some practical details about how to use PCV are given below.

## How to use PCV in R

We created a dedicated package, `pcv`, which contains all PCV functions. You can download the package file from [Releases](https://github.com/svkucheryavski/pcv/releases) and install from the file. The package is currently under review in CRAN and as soon as it is accepted you will be also able install the package from CRAN by running `install.packages("pcv")`.

There are four main functions in the package:

* `pcv()` is implementation of the original version, we keep it for compatibility, if you consider using PCV a new project, do not use it, use the new version.

* `pcvpca()` is implementation of the new version for PCA/SIMCA models.

* `pcvpcr()` is implementation of the new version for PCR models.

* `pcvpls()` is implementation of the new version for PLS models.

All four functions return PV-set generated with a given parameters, which has the same size as the calibration set. In case of regression PV-set is generated only for predictors (X), the response values for PV-set are the same as for the calibration set.

The last two functions return the PV-set with additional attribute, `"D"` which is matrix containing scaling factors ($c_k/c$), for each segment and each component. See all details in the paper [3]. The matrix can be visualized as a heatmap, similar to shown in the paper, using method `plotD()` which is also a part of the package.

Below are examples of the functions syntax with all parameters:

```r
# for PCA/SIMCA models
Xpv <- pcvpca(X, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PCR models
Xpv <- pcvpcr(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PLS models
Xpv <- pcvpls(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# get matrix D and show its values as boxplot
D <- attr(Xpv, "D")
boxplot(D)

# show heatmap for D values
plotD(Xpv)
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix, not a data frame). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. As mentioned above, the method generates PV-set only for predictors, the response values for the calibration set and for the PV-set are the same.

Parameter `ncomp` is a number of principal components in case of PCA/PCR models or number of latent variables in case of PLS based method. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `center` and `scale` define if the predictors must be mean centered and/or standardized. By default `center = TRUE` and `scale = FALSE`. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with the calibration set.

Finally, parameter `cv` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `cv = 4`). In this case this number specifies number of segments for random splits, except `cv = 1`, which is a special case for leave-one-out (full cross-validation).

* A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in the code examples above, tells PCV to use Venetian blinds splits with 4 segments.

* A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number of values in this vector must be the same as number of rows in the training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.


File `demo.R` from this repository contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.


## How to use in MATLAB

All functions for running PCV in MATLAB are combined into a dedicated Toolbox which is available from Mathworks File Exchange repository ([direct link](https://se.mathworks.com/matlabcentral/fileexchange/121468-procrustes-cross-validation)). Alternatively you can click on "Add-Ons" button in MATLAB, select "Get Add-Ons", and search "Procrustes cross-validation".

You can also download the toolbox as `.mtlbx` file from [Releases](https://github.com/svkucheryavski/pcv/releases) section and install it from file.

The Toolbox has Getting started tutorial, where you can see the examples of code and try them.

Below you will also find a very short description of main functions:

```matlab
% for PCA/SIMCA models
Xpv = pcvpca(X, 20, true, false, {"ven", 4});

% for PCR models
Xpv = pcvpcr(X, y, 20, true, false, {"ven", 4});

% for PLS models
Xpv = pcvpls(X, y, 20, true, false, {"ven", 4});
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. The method generates PV-set only for predictors, the response values for calibration set and PV-set are the same.

The next parameter (shown as 20 in the examples above) is a number of principal components in case of PCA/PCR models or number of latent variables in PLS model. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Next two parameters with logical values define if the predictors must be mean centered and/or standardized. By default `Center` is set to `true` and `Scale` is set to  `false` as in the examples above. Regardless which settings you use the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with calibration set.

Finally, the last parameter defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `4`). In this case this number specifies number of segments for random splits, except `1` which is a special case for leave-one-out (full cross-validation).

* A cell array with 2 values: `{"name", nseg}`. In this case `"name"` defines the way to make the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into.

* A vector with integer numbers, e.g. `[1, 2, 3, 1, 2, 3, 1, 2, 3]`. In this case number of values in this vector must be the same as number of rows in training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.

File `demo.m` contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.

## References

Papers [1] and [2] describe the original version, paper [3] — the new version.

1. Kucheryavskiy, S., Zhilin, S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation—A Bridge between Cross-Validation and Independent Validation Sets.* Analytical Chemistry,  92 (17), 2020. pp.11842–11850. DOI: [10.1021/acs.analchem.0c02175](https://doi.org/10.1021/acs.analchem.0c02175)

2. Pomerantsev A. L., Rodionova O. Ye. *Procrustes Cross-Validation of short datasets in PCA context.* Talanta, 226, 2021. DOI: [10.1016/j.talanta.2021.122104](https://doi.org/10.1016/j.talanta.2021.122104)

3. Kucheryavskiy S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation of multivariate regression model*. Submitted, 2022.


## Bugs and improvements

The code will be improved and extended gradually. If you found a bug please report using issues or send an email.


