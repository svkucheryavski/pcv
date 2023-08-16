# Procrustes cross-validation

<img src="logo.png" width="75" height="75" style="float:left;padding-right:10px;"> Procrustes cross-validation (PCV) is a novel approach for validation of chemometric models. It generates a new dataset — *PV-set*, which can be used for validation of PCA/SIMCA/PCR/PLS models in the same way as with an independent validation set.

You can learn more about the method in [this paper](https://doi.org/10.1016/j.aca.2023.341096). The paper is Open Access and is freely available for everyone. For more information see *References* section below.

This repository contains source code of PCV implementation in several programming languages (R, MATLAB, Python and Javascript) as well as some practical details. If you do not use any of the languages you can apply PCV to your own data via an interactive web-application: [https://mda.tools/pcv/](https://mda.tools/pcv/). The application takes your data as a CSV file (in case of regression first column must contain the response values), it does all calculations in your browser and does not send any information to server.

## Release notes

Last minor update [1.1.x) was released on 12th of August, 2023 for all languages. It contains small improvements, better test coverage and new experimental feature — CV-scope. This feature lets you define how local calibration sets must be centered and scaled.

By default a *global* scope is used, so all local models have the same center as the global model in the original variable space. In this case all local calibration sets are autoscaled using globally computed mean and standard deviation values. This implementation is also described in [the paper](https://doi.org/10.1016/j.aca.2023.341096).

If you use a local CV-scope, then each local calibration and validation sets will be autoscaled independently by using mean values and standard deviations computed for the local calibration set. This will slightly increase the sampling error and, in theory, should make it more realistic. But the Procrustean rules for PV-sets computed using the local CV scope will only hold approximately true (there is a strong agreement in case of the global scope). This feature requires additional investigation so use it with caution.

Finally I have also implemented PCV for Python and Javascript languages, check links with additional information below.


## Implementation of PCV in different languages

Please use the links below to get to a part of the repository, describing implementation of PCV in a particular language:

* [R](R/README.md)
* [MATLAB](MATLAB/README.md)
* [Python](Python/README.md)
* [Javascript](Javascript/README.md)


## References

Here you can find some useful papers about PCV. Paper [1] describes the initial version of the method. Paper [2] shows some practical examples for using PCV in SIMCA models. Paper [3] introduces the new, more general and efficient version of the method.

Please cite the paper [3] if you use PCV in your projects.

1. Kucheryavskiy, S., Zhilin, S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation—A Bridge between Cross-Validation and Independent Validation Sets.* Analytical Chemistry,  92 (17), 2020. pp.11842–11850. DOI: [10.1021/acs.analchem.0c02175](https://doi.org/10.1021/acs.analchem.0c02175)

2. Pomerantsev A. L., Rodionova O. Ye. *Procrustes Cross-Validation of short datasets in PCA context.* Talanta, 226, 2021. DOI: [10.1016/j.talanta.2021.122104](https://doi.org/10.1016/j.talanta.2021.122104)

3. Kucheryavskiy S., Rodionova, O., & Pomerantsev, A. *Procrustes Cross-Validation of multivariate regression model*. Analytica Chimica Acta, 1255, 2023 [10.1016/j.aca.2023.341096](https://doi.org/10.1016/j.aca.2023.341096).


## Bugs and improvements

The method and the implementations will be improved and extended gradually. If you found a bug please report using [issues](https://github.com/svkucheryavski/pcv/issues) or send an email.


