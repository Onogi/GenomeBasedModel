# GenomeBasedModel
The package has been developed for estimating parameters of mathematical biological models and effects of genome-wide SNPs on the model parameters simultaneously by building prior distributions based on the SNPs for the parameters.

The estimation framework is based on a variational Bayes-MCMC hybrid algorithm. The prototype of the framework can be found in Onogi et al. (2016) Theor. Appl. Genet. 129:805-817 (https://link.springer.com/article/10.1007%2Fs00122-016-2667-5). The generalized algorithm and its superiority over existing methods are illustrated in Onogi (2020) Bioinformatics (https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa129/5758262). The paper of Bioinformatics is not open access, but it's possible to share privately through researchgate (https://www.researchgate.net/profile/Akio_Onogi). Please contact through researchgate if you have interests.

The latest version is 1.8. Several bugs were fixed, and imputation of missing observations were quited.

For installation, Rcpp and rrBLUP packages are required (the latter is for creating the genetic relationship matrix). Rcpp and rrBLUP versions of the tested environment are 1.0.7 and 4.6, respectively.

Usages of the package are illustrated in the vignette of the package. Type ?GenomeBasedModel after loading the package. Examples also can be found in https://github.com/Onogi/IntegratingCGMwithGP which includes scripts used in "Integration of crop growth models and genomic prediction" by Akio Onogi, a chapter of "Genomic prediction of complex traits" editted by Nour Ahmadi and Jerome Bartholomé. Here GenomeBasedModel is applied to two crop growth models (DVR model and maize growth model). R scripts in AnalysisExample_DVRmodel and AnalysisExample_MaizeGrowthModel will be heplful.

Scripts_Onogi(2020)Bioinformatics.zip includes all R and Cpp scripts rquired for the experiments in Onogi (2020) Bioinformatics which also will be good illustrations of the usage of the package.

Archive includes previous versions.
