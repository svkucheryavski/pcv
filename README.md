# Procrustes cross-validation

This repository contains implementation of Procrustes cross-validation in R and MATLAB. The current version of code is *0.1.1*. The code can be also downloaded as archive from [Release](https://github.com/svkucheryavski/pcv/releases) section.

## Short description

Procrustes cross-validation is a new approach for validation of chemometric models. It makes possible to generate a new dataset, named *pseudo-validation set* and use it for validation of models in the same way as with an independent test set. The generation is done using following algorithm:

1. A global PCA model is created using calibration set **X** and *A* components
2. The rows of matrix **X** are split into *K* segments using venetian blinds splitting
3. For each segment *k* from *{1, 2, ..., K}*:
    * a local PCA model is created using the rows from the segment
    * an angle between the local and the global model is estimated
    * rows from the current segment are rotated in original variable space by the angle
4. All rotated measurements are then combined into a matrix with pseudo-validation set

So, pseudo-validation set is built on top of the calibration set but has its own sampling error. Since it is not independent from the calibration set we recommend to limit its use by model optimization and do not use it for assessment of performance of the final model. More details can be found in our [preprint](), the paper is currently under review in one of the scientific journals.

## How to use in R

All code is located in `pcv.R`, no extra packages are needed. Simply copy the file to your working directory and use `source("pcv.R")` inside your script to load all necessary functions to the environment. Then the syntax is follows:

```r
pvset <- pcv(X, ncomp = 10, nseg = 4, scale = FALSE)
```

Here `X` is your calibration set (as a numerical matrix), `ncomp` is a number of prinicpal components for PCA decomposition (use at lease enough to explain 95-99% of variation of the data values), `nseg` is the number of segments for cross-validation procedure. So far only systematic cross-validation (venetian blinds) is supported, so make sure that rows of `X` are sorted correctly or shuffled. Parameter `scale` allows to standardize your data prior to the generation, which is useful if your variables have different nature and/or units. The generated data will be unscaled though.

File `demo.R` contains a demo code based on $NIRSim$ dataset from the paper. See comments in the code for more details.

## How to use in MATLAB

See details for R implementation above.  To use it in MATLAB you have to copy file `pcv.m` to your script folder (or any other folder MATLAB knows path to). Then use the following syntax:

```matlab
pvset = pcv(X, 'nComp', 4, 'nSeg', 4, 'Scale', true);
```

File `demo.m` contains a demo code based on $NIRSim$ dataset from the paper. See comments in the code for more details.

### Bugs and improvements

We plan to improve and extend the code gradually. If you found a bug please report using issues or send an email.




