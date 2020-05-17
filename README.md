# Procrustes cross-validation

This repository contains implementation of Procrustes cross-validation in R and MATLAB. The current version of code is 0.1.1. The code can be also downloaded as archive from [Release]() section. 

## Short description

Procrustes cross-validation is a new approach for validation of chemometric models. It allows you to generate a new dataset, named *pseudo-validation set* and use it for validation of your models in the same way one usually uses an independent test or validation set. The generation is done using following algorithm:

1. A global PCA model is created using calibration set X and A components
2. The rows of matrix X are split into K segments using venetian blinds spliting
3. For each segment k:
    * a local PCA model is created using the rows from the segment
    * an angle between the local and the global model is estimated
    * rows from the current segment are rotated in original variable space by the angle
4. All rotated measurements are then combined into a matrix with pseudo-validation set

So, pseudo-validation set is built on top of the calibration set but has its own sampling error. More details can be found in 

## How to use in R

All code is located in `pcv.R`, no extra packages are needed. Simply copy the file to your working directory and use `source("pcv.R")` inside your script to load all necessary functions. Then the syntax is follows:

```r
pvset <- pcv(X, A = 4, K = 4, scale = TRUE)
```

Here `X` is your calibration set (as a numerical matrix), `A` is a number of prinicpal components for PCA decomposition (use at lease enough to explain 95-99% of variation of the data values), `K` is number of segments for cross-validation procedure. So far only systematic cross-validation (venetian blinds) is supported, so make sure that rows of `X` are sorted correctly or shuffled. Parameter `scale` allows to standardize your data prior to generation, which is usefull if your variables have different nature and/or units.

## How to use in MATLAB

See details for R implementation above.  To use it in MATLAB you have to copy file `pcv.m` to your script folder (or any other folder MATLAB knows path to). Then use the following syntax:

```matlab
pvset = pcv(X, 'nComp', 4, 'nSegments', 4, 'Scale', true);
```



