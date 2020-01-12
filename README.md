# GenomeBasedModel
This is the page for an R package GenomeBasedModel.
The package has been developed for estimating parameters of mathematical biological models and effects of genome-wide SNPs on the model parameters simultaneously by building prior distributions based on the SNPs for the parameters.

The estimation framework is based on a variational Bayes-MCMC hybrid algorithm.
The prototype of the framework can be found in Onogi et al. (2016) Theor. Appl. Genet. 129:805-817.
The package and its performance are introduced in Onogi (2020) Bioinformatics (under review).

The latest version is 1.3.

tar.gz files are the source packages and .zip are the binary packages.

Instration can be done following the standard procedure for R package instration.
Rcpp and rrBLUP packages are required (the latter is for creating the genetic relationship matrix).
Rcpp and rrBLUP versions of the tested environment are 1.0.0 and 4.6, respectively.

Usage of the package is illustrated in the R documentations of the package.
Please see the examples in GenomeBasedModel function of the package for the usage.
