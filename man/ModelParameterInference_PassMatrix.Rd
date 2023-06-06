\name{ModelParameterInference_PassMatrix}
\alias{ModelParameterInference_PassMatrix}
\title{Inferring model parameters using MCMC}
\description{
This function is called by GenomeBasedModel when PassMatrix is TRUE for inferring model parameters using MCMC. Although it is possible to use this function to infer parameters of arbitrary models, no warning is given for irregular arguments because this C++ function is designed for internal use.
}
\usage{
ModelParameterInference_PassMatrix(Input, Freevec, Y, Missing, Y_uncertainty,
                                   Nl, Ne, Np, Nr, SampledPara, SampledVe, MeanOfPrior,
                                   VarOfPrior, Ve, Nve, Residualgroup, SdforParameters,
                                   SdforVe, Loglike, Burnin, Thi, Totalite, Order,
                                   Lowerlimit, Upperlimit, AccParameters, AccVe, Model,
                                   NeNl, Transformation, nRandom_para, nRandom_ve,
                                   uRandom_para, uRandom_ve)
}
\arguments{
\item{Input}{Same as Input of GenomeBasedModel.}
\item{Freevec}{Same as Freevec of GenomeBasedModel}
\item{Y}{An Ne x Nl matrix of observations where Ne is the dimension of observations and Nl is the number of genotypes. Genotype IDs, which are required by GenomeBasedModel, are not required here. The likelihood for genotype i is calculated from Y[,i] and Modelfunction(Input[,i], Freevec, parameters).}
\item{Missing}{A scalar to specify missing values in Y (NA IS NOT ALLOWED).}
\item{Y_uncertainty}{An Ne x Nl matrix of uncertainties of Y. These values are added to the residual variance.}
\item{Nl}{The number of genotypes (the number of columns of Y).}
\item{Ne}{The number of rows of Y}
\item{Np}{The number of parameters in Model function}
\item{Nr}{The number of rows of Input}
\item{SampledPara}{An (Np x (Totalite-Burnin)/Thi) x Nl matrix for saving MCMC samples of model parameters. Fill the elements with arbitrary numbers, for example, 0.}
\item{SampledVe}{A (Totalite-Burnin)/Thi x Nve matrix for saving MCMC samples of residual variances. Fill the elements with arbitrary numbers, for example, 0.}
\item{MeanOfPrior}{An Np x Nl matrix specifying the means of prior distributions of model parameters. This implies that this function assumes independent prior distributions among genotypes. In GenomeBasedModel, however, dependence arises because the prior distributions are connected using genome-wide SNPs.}
\item{VarOfPrior}{An Np x Nl matrix specifying the variances of prior distributions of model parameters.}
\item{Ve}{An Nve length vector specifying the residual variances. This vector is updated during optimization. Fill the elements with numbers > 0 (e.g., 1).}
\item{Nve}{The number of residual variance groups}
\item{Residualgroup}{An Ne length vector specifying the residual variance group. Note that the index should start from 0.}
\item{SdforParameters}{Same as SdforParameters in GenomeBasedModel}
\item{SdforVe}{Same as SdforVe in GenomeBasedModel}
\item{Loglike}{A Totalite/Thi length vector for saving log likelihood. Fill the elements with arbitrary numbers, for example, 0.}
\item{Burnin}{A scalar specifying the number of iterations discarded as burnin.}
\item{Thi}{A scalar specifying the sampling interval.}
\item{Totalite}{A scalar specifying the number of MCMC iterations. The number of MCMC samples is Totalite/Thi.}
\item{Order}{An Nl length vector specifying the update order of genotypes. Different orders may show different convergence, but this argument is fixed to 0:(Nl-1) in GenomeBasedModel. Note that the index should start from 0.}
\item{Lowerlimit}{An Np length vector to specify the lower limits of parameters.}
\item{Upperlimit}{An Np length vector to specify the upper limits of parameters.}
\item{AccParameters}{An Np length vector with elements being 0 [i.e., rep(0, NP)]. This is for counting the acceptance of parameter updates in MCMC.}
\item{AccVe}{An Nve length vector with elements being 0 [i.e., rep(0, Nve)]. This is for counting the acceptance of parameter updates in MCMC.}
\item{Model}{Same as Modelfunction in GenomeBasedModel.}
\item{NeNl}{An Nve length vector specifying the number of genotypes with non-missing records in Y for each residual group. This is used for update of the residual variances}
\item{Transformation}{An Np length vector specifying parameter transformation methods in Modelfunction 0:not done, 1:log, and 2:logit}
\item{nRandom_para}{A (Totalite + 1) x (Np x Nl) matrix of random variables from N(0, 1). Used for MCMC update of model parameters.}
\item{nRandom_ve}{A Totalite x Nve matrix of random variables from N(0, 1). Used for MCMC update of residual variance.}
\item{uRandom_para}{A Totalite x (Np x Nl) matrix of random variables from Uniform(0, 1). Used for MCMC update of model parameters.}
\item{uRandom_ve}{A Totalite x Nve matrix of random variables from Uniform(0, 1). Used for MCMC update of residual variance.}
}
\value{
Return 0 when finished successfully. The results of MCMC are written in SampledPara and SampledVe. In SampledPara, MCMC samples for each parameter are vertically stacked. That is, SampledPara[1,1] includes the first MCMC sample of the first parameter for the first genotype and SampledPara[1+(Totalite-Burnin)/Thi,1] includes the first MCMC sample of the second parameter for the first genotype.
}
\references{
  Onogi A. Connecting mathematical models to genomes: Joint estimation of model parameters and genome-wide marker effects on these parameters. Bioinformatics 36, 3169-3176 (2020).\cr
}
\author{
  Akio Onogi
}

