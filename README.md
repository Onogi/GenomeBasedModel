# GenomeBasedModel
The package has been developed for estimating parameters of mathematical biological models and effects of genome-wide SNPs on the model parameters simultaneously by building prior distributions based on the SNPs for the parameters.

The estimation framework is based on a variational Bayes-MCMC hybrid algorithm.
The prototype of the framework can be found in Onogi et al. (2016) Theor. Appl. Genet. 129:805-817.
The generalized algorithm and its superiority over existing methods are illustrated in Onogi (2020) Bioinformatics (https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa129/5758262).
The paper of Bioinformatics is not open access, but it's possible to share the paper privately through researchgate (https://www.researchgate.net/profile/Akio_Onogi).
Please contact through researchgate if you have interests.

The latest version is 1.3.
tar.gz files are the source packages and .zip are the binary packages.

Instration can be done following the standard procedure for R package instration.
Rcpp and rrBLUP packages are required (the latter is for creating the genetic relationship matrix).
Rcpp and rrBLUP versions of the tested environment are 1.0.0 and 4.6, respectively.

Usage of the package is illustrated in the R documentations of the package.
Please see the examples in GenomeBasedModel function of the package for the usage.

Scripts_Onogi(2020)Bioinformatics.zip includes all R and Cpp scripts rquired for the experiments in Onogi (2020) Bioinformatics which also will be good illustrations of the usage of the package.
