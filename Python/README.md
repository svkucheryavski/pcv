# Procrustes cross-validation for Python

Package `pcv` implements [Procrustes cross-validation](https://github.com/svkucheryavski/pcv) in Python language.

Last version of the package (*1.0.0*) was released 8th of August, 2023 and contains small improvements, better test coverage, as well as a new experimental feature — CV scope. See details in the overall [project description](https://github.com/svkucheryavski/pcv).

## Getting started

You can install the package directly from PyPI by running `pip install pcv`. If you already have the package installed and want to update it to the newest version use: `pip install pcv --upgrade`

There are three main functions in the package:

* `pcvpca()` is implementation of PCV for PCA/SIMCA models.

* `pcvpcr()` is implementation of PCV for PCR models.

* `pcvpls()` is implementation of PCV for PLS models.

All three functions return PV-set generated with given parameters. The PV-set has the same size as the calibration set. In case of regression (PCR or PLS) PV-set is generated only for predictors (X), the response values for PV-set are the same as for the calibration set.

The last two functions return the PV-set with additional attribute, `"D"` which is matrix containing scaling factors ($c_k/c$), for each segment and each component. See all details in the [paper](https://doi.org/10.1016/j.aca.2023.341096). The matrix can be visualized as a heatmap, similar to the ones shown in the paper, using method `plotD()` which is also a part of the package.

Below are examples of the functions syntax with all parameters:

```python
# for PCA/SIMCA models
Xpv = pcvpca(X, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PCR models
Xpv = pcvpcr(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# for PLS models
Xpv = pcvpls(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))

# show heatmap for D values
plotD(Xpv)

# get the matrix D and show its values as boxplot
D <- attr(Xpv, "D")
boxplot(D)
```

Here `X` is a matrix with predictors for your calibration set (as a numerical matrix, not a data frame). In case of regression model you also need to provide a vector or a matrix with response values for the training set, `y`. As mentioned above, the method generates PV-set only for predictors, the response values for the calibration set and for the PV-set are the same.

Parameter `ncomp` is a number of principal components in case of PCA/PCR models or number of latent variables in case of PLS based method. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `center` and `scale` define if the predictors must be mean centered and/or standardized. By default `center = TRUE` and `scale = FALSE`. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with the calibration set.

Finally, parameter `cv` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `cv = 4`). In this case this number specifies number of segments for random splits, except `cv = 1`, which is a special case for leave-one-out (full cross-validation).

* A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in the code examples above, tells PCV to use Venetian blinds splits with 4 segments.

* A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number of values in this vector must be the same as number of rows in the training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.

As it is written above, from *1.1.0*, there is additional parameter, `cv.scope`, which can have one of the two values, `"global"` or `"local"`. The default value is `"global"`, if you want to try the local scope, just add this parameter when you call one of the functions, like shown below:


```r
# PCV for PLS models with local CV scope
Xpv <- pcvpls(X, y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4), cv.scope = "local")
```

File `demo.R`, which you can download from this repository contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.

The package code will be improved and extended gradually. If you found a bug please report using [issues](https://github.com/svkucheryavski/pcv/issues) or send an email.


