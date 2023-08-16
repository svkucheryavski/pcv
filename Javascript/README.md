# Procrustes cross-validation for Javascript

Package `pcv` implements [Procrustes cross-validation](https://github.com/svkucheryavski/pcv) in Javascript language.

Last main version of the package was released in August, 2023 and contains small improvements, better test coverage, as well as a new experimental feature — CV scope. See details in the overall [project description](https://github.com/svkucheryavski/pcv).

## Getting started

You can install the package directly from NPM by running `npm install pcv`. There are three main functions in the package:

* `pcvpca()` is implementation of PCV for PCA/SIMCA models.

* `pcvpcr()` is implementation of PCV for PCR models.

* `pcvpls()` is implementation of PCV for PLS models.

All three functions return PV-set generated with given parameters. The PV-set has the same size as the calibration set. In case of regression (PCR or PLS) PV-set is generated only for predictors (X), the response values for PV-set are the same as for the calibration set.

All outcomes are instances of class `Matrix` from [mdatools](https://www.npmjs.com/package/mdatools) package, check its description for more details. It also uses implementation of PCA, PCR and PLS from the package.

The last two functions return the PV-set and an additional outcome, `D`, which is a matrix containing scaling factors ($c_k/c$), for each segment and each component. See all details in [the paper](https://doi.org/10.1016/j.aca.2023.341096).

Below are examples of the function syntax with main parameters:

```javascript
import { pcvpca, pcvpcr, pcvpls } from 'pcv';
import { pcafit, pcapredict, pcrfit, pcrpredict, plsfit, plspredict } from 'mdatools/models';

// set cross-validation settings — common for all methods above
const cv = {type: 'ven', nseg: 4};

// for PCA/SIMCA models:
// create a global PCA model with A = 20, centring but not scaling
// then we crate PV-set based on the model
const mpca = pcafit(X, 20, true, false);
const Xpv_pca = pcvpca(X, mpca, 20, cv);

// for PCR models:
// create a global PCR model with A = 20, centring but not scaling
// then we get a PV-set and matrix with scalars based on the model
const mpcr = pcrfit(X, Y, 20, true, false);
const [Xpv_pcr, D_pcr] = pcvpcr(X, Y, mpcr, 20, cv);

// for PLS models:
// create a global PLS model with A = 20, centring but not scaling
// then we get a PV-set and matrix with scalars based on the model
const mpls = plsfit(X, Y, 20, true, false);
const [Xpv_pls, D_pls] = pcvpls(X, Y, mpls, 20, cv);
```

Here `X` is a matrix with predictors for your calibration set (instance of `Matrix` class from `mdatools` package). In case of regression model you also need to provide a matrix with response values for the training set, `Y`. The method generates PV-set only for predictors, the response values for the calibration set and for the PV-set are the same.

Parameter `ncomp` is a number of principal components in case of PCA/PCR models or number of latent variables in case of PLS based method. Number of components must be large enough, larger than the expected optimal number. In case of PCA use components which explain at least 99% of the data.

Parameters `center` and `scale` define if columns of `X` and `Y` should be mean centered and/or standardized. By default `center = true` and `scale = false`. Regardless which settings you use, the resulted PV-set will be in original units (uncentered and unstandardized), so you can compare it directly with the calibration set.

Finally, parameter `cv` defines how to split the rows of the training set. The split is similar to cross-validation splits, as PCV is based on cross-validation resampling. This parameter can have the following values:

* A JSON with 2 field: `'type'` and `'nseg'`. In this case `'type'` defines the way to make the splits. You can select one of the following: `'loo'` for leave-one-out, `'rand'` for random splits or `'ven'` for Venetian blinds (systematic) splits. The second field, `'nseg'`, should be a number of segments for splitting the rows into. For example, `cv = {type: 'ven', nseg: 4}`, shown in the code examples above, tells PCV to use Venetian blinds splits with 4 segments.

* A vector with integer numbers, e.g. `cv = index([1, 2, 3, 1, 2, 3, 1, 2, 3])`. Here `index` is a function from `mdatools` package. In this case number of values in this vector must be the same as number of rows in the training set. The values specify which segment a particular row will belong to. For the example shown here, it is assumed that you have 9 rows in the calibration set, which will be split into 3 segments. The first segment will consist of measurements from rows 1, 4 and 7.

There is also additional parameter, `cvscope`, which can have one of the two values, `'global'` or `'local'`. The default value is `'global'`, if you want to try the local scope, just add this parameter when you call one of the functions, like shown below (in this example some of the arguments are skipped for simplicity):

```javascript
// PCV for PLS models with local CV scope
const [Xpv_pls, D_pls] = pcvpls(X, Y, mpls, 20, cv, 'global');
```

The package code will be improved and extended gradually. If you found a bug please report using [issues](https://github.com/svkucheryavski/pcv/issues) or send an email.


