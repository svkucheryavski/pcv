# Procrustes cross-validation for R

Package `pcv` implements [Procrustes cross-validation](https://github.com/svkucheryavski/pcv) in R language.

Last version of the package (*1.2.0*) was released in December, 2025 and contains small improvements to avoid issues with NaN numbers in case if original dataset or local subsets contain columns with constant values. See details in the overall [project description](https://github.com/svkucheryavski/pcv).

## Getting started

You can install the package directly from CRAN by running `install.packages("pcv")`. If you already have the package installed and want to update it to the newest version use: `update.packages("pcv")`

There three main functions in the package:

* `pcvpca()` is implementation of PCV for PCA/SIMCA models.

* `pcvpcr()` is implementation of PCV for PCR models.

* `pcvpls()` is implementation of PCV for PLS models.

All three functions return PV-set generated with given parameters. The PV-set has the same size as the calibration set. In case of regression (PCR or PLS) PV-set is generated only for predictors (X), the response values for PV-set are the same as for the calibration set.

The last two functions return the PV-set with additional attribute, `"D"` which is matrix containing scaling factors ($c_k/c$), for each segment and each component. See all details in the [paper](https://doi.org/10.1016/j.aca.2023.341096). The matrix can be visualized as a heatmap, similar to the ones shown in the paper, using method `plotD()` which is also a part of the package.

Below are examples of the functions syntax with all parameters:

```r
# for PCA/SIMCA models
Xpv <- pcvpca(X, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PCR models
Xpv <- pcvpcr(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PLS models
Xpv <- pcvpls(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# show heatmap for D values
plotD(Xpv)

# get the matrix D and show its values as boxplot
D <- attr(Xpv, "D")
boxplot(D)
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix, not a data frame). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. As mentioned above, the method generates PV-set only for predictors, the response values for the calibration set and for the PV-set are the same.

Parameter `ncomp` is a number of principal components in case of PCA/PCR models or number of latent variables in case of PLS based method. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `center` and `scale` define if the predictors must be mean centered and/or standardized. By default `center = TRUE` and `scale = FALSE`. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with the calibration set.

Finally, parameter `cv` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation resampling. This parameter can have the following values:

* A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in the code examples above, tells PCV to use Venetian blinds splits with 4 segments.

* A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number of values in this vector must be the same as number of rows in the training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.

As it is written above, from *1.1.0*, there is additional parameter, `cv.scope`, which can have one of the two values, `"global"` or `"local"`. The default value is `"global"`, if you want to try the local scope, just add this parameter when you call one of the functions, like shown below:


```r
# PCV for PLS models with local CV scope
Xpv <- pcvpls(X, y, ncomp = 20, cv = list("ven", 4), cv.scope = "local")
```

File `demo.R`, which you can download from this repository contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.

The package code will be improved and extended gradually. If you found a bug please report using [issues](https://github.com/svkucheryavski/pcv/issues) or send an email.


## Three-way data

From version *1.3.0* the package also implements PCV for three-way data, for example fluorescence excitation-emission matrices (EEM), where every sample is a matrix (e.g. excitation wavelengths x emission wavelengths). The calibration set is a three-way array with samples in the first mode. There are three methods:

* `pcv3wfu()` unfolds every sample to a vector and generates PV-set with PCA based PCV (full unfolding, FU).

* `pcv3wbu()` splits the array into blocks along mode 2 or mode 3 and generates PV-set for each block with PCA based PCV, using the same split of the samples (block unfolding, BU).

* `pcv3whsv()` is based on HO-SVD (Tucker2) model with components in modes 2 and 3 (HSV). Projected onto the global HO-SVD model, the PV-set gives the same scores and orthogonal distances as cross-validation of this model.

All three functions return PV-set as a three-way array of the same size as the calibration set, in the original units:

```r
# full unfolding
Xpv <- pcv3wfu(X, ncomp = 20, cv = list("ven", 4))

# block unfolding along mode 2
Xpv <- pcv3wbu(X, mode = 2, ncomp = 10, cv = list("ven", 4))

# HO-SVD based method with 6 x 6 components
Xpv <- pcv3whsv(X, ncomp = c(6, 6), cv = list("ven", 4))
```

Parameter `cv` works the same way as for the two-way methods. All three functions also have parameter `seed`, which makes random splits and random parts of the PV-set reproducible.

Two practical recommendations:

* HSV PV-set is a good choice in general. For validation of PARAFAC and Tucker3 models it can be replaced with FU PV-set.

* Generate PV-set with more components than the largest model you are going to validate. A PV-set carries the sampling variation only within the subspace it was generated in.

File `demo3w.R` contains a demo code with simulated EEM data, which shows how to generate the PV-sets and use them for validation of a three-way model. The methods are described in: S. Kucheryavskiy, *Procrustes cross-validation for three-way data* (submitted).
