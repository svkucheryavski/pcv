# Procrustes cross-validation for MATLAB

This part of repository contains implementation of [Procrustes cross-validation](https://github.com/svkucheryavski/pcv) in MATLAB as a dedicated toolbox.

Last version of the toolbox (*1.1.0*) was released 8th of August, 2023 and contains small improvements, better test coverage, as well as a new experimental feature — CV scope. See details about this feature in the overall [project description](https://github.com/svkucheryavski/pcv).

## Getting started

All functions for running PCV in MATLAB are combined into a dedicated Toolbox which is available from Mathworks File Exchange repository ([direct link](https://se.mathworks.com/matlabcentral/fileexchange/121468-procrustes-cross-validation)). Alternatively you can click on "Add-Ons" button in MATLAB, select "Get Add-Ons", and search "Procrustes cross-validation".

You can also download the toolbox as `.mtlbx` file from [Releases](https://github.com/svkucheryavski/pcv/releases) section and install it from file.

The Toolbox has "Getting started" tutorial, where you can see the examples of code and try them.

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

Next two parameters with logical values define if the predictors must be mean centered and/or standardized. By default `Center` is set to `true` and `Scale` is set to  `false` as in the examples above. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with calibration set.

Finally, the last parameter defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation. This parameter can have the following values:

* A number (e.g. `4`). In this case this number specifies number of segments for random splits, except `1` which is a special case for leave-one-out (full cross-validation).

* A cell array with 2 values: `{"name", nseg}`. In this case `"name"` defines the way to make the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a number of segments for splitting the rows into.

* A vector with integer numbers, e.g. `[1, 2, 3, 1, 2, 3, 1, 2, 3]`. In this case number of values in this vector must be the same as number of rows in training set. The values specify which segment a particular row will belong to. In case of the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.

As it is written above, from *1.1.0*, there is an additional parameter, `CVScope`, which can have one of the two values, `"global"` or `"local"`. The default value is `"global"`, if you want to try the local scope, just add this parameter when you call one of the functions, like shown below:


```matlab
% PCV for PLS models with local CV scope
Xpv = pcvpls(X, y, 20, true, false, {"ven", 4}, "local");
```

File `demo.m`, which you can download from this repository contains a demo code based on *Corn* dataset from the paper to be published. See comments in the code for more details.

The code will be improved and extended gradually. If you found a bug please report using [issues](https://github.com/svkucheryavski/pcv/issues) or send an email.


