\name{GenomeBasedModel}
\alias{GenomeBasedModel}
\title{Building genome-based prior distributions on statistical models in biology}
\description{
This function builds prior distributions for the parameters of statistical models in biology that are assumed to be genotype-specific (individual-specific) using genome-wide SNPs. The parameters of the model are regressed on genome-wide SNPs and the parameters and the SNP effects on the parameters are inferred jointly. Specifically, the function repeats model parameter inference (MPI) and genome-wide regression (GWR) for SNP effect estimation until convergence. The built genome-based model is expected to avoid overfitting and to be robust and accurate because it uses observations from all genotypes to estimate genotype-specific parameters.
}
\usage{
GenomeBasedModel(Input, Freevec, Y, Missing, Np, Geno, Methodcode, Referencevalues,
                Modelfunction, PolygenicEffect = FALSE,
                Transformation = NULL, Residualgroup = NULL, Hyperpara = NULL,
                Q = NULL, SdforParameters = NULL, SdforVe = NULL,
                K = NULL, Y_uncertainty = NULL, Nite_GWR = 1000, Nloop = 200,
                Nite_MPI = 500, Nite_MPI_last = 3000, Thin=10,
                Tol = 1e-9, Nblock = 20, Rhatthr = 1.01, StopN = 10,
                Lowerlimit = rep(-1e+20, Np), Upperlimit = rep(1e+20, Np),
                Initialvalues = NULL, PassMatrix = FALSE)
}
\arguments{
  \item{Input}{
  A matrix of which the number of columns is the number of genotypes (Ny). The number of rows is arbitrary. Input is passed to Modelfunction, which describes the model that the user wants to build a genome-based model.
}
  \item{Freevec}{
  An arbitrary length vector of integers. This vector is passed to Modelfunction and can be used for any purpose, for example, to define partitions in Input to put multiple information in Input. ENTER AN ARBITRARY INTEGER (e.g., 0) IF NOT USED.
}
  \item{Y}{
  A (D + 1) x Ny matrix of observations where D is the dimension of observations. THE FIRST LINE INCLUDES THE GENOTYPE IDs (only integers). Missing values in Y are specified by Missing. Likelihood is calculated from Y[-1,] and outputs of Modelfunction.
}
  \item{Missing}{
  A scalar to specify missing values in Y (NA IS NOT ALLOWED).
}
  \item{Np}{
  An integer to specify the number of parameters used in Modelfunction.
}
  \item{Geno}{
  A (P + 1) x Nx matrix of SNP genotypes where P is the number of SNPs and Nx is the number of genotypes with SNP genotypes. THE FIRST LINE INCLUDES THE GENOTYPE IDs (only integers). Genotypes not included in Input and Y can be included (thus, Nx can be > Ny). Note that the function does nothing for Geno. Thus, standardization or genotype imputation should be performed by the user, if necessary. It is also notable that, when the user uses GBLUP and K is NULL, SNP genotypes should be coded as -1 (aa), 0 (Aa), and 1 (AA), which are required by A.mat of rrBLUP called for genetic relationship matrix computation.
}
  \item{Methodcode}{
  An Np length vector to specify regression methods in GWR for each parameter. The regression methods are specified with codes:
    \itemize{
    \item 1: Bayesian Lasso (BL)
    \item 2: extended Bayesian Lasso (EBL)
    \item 3: weighted Bayeisan shrinkage regression (wBSR)
    \item 4: BayesC
    \item 5: Stochastic search variable selection (SSVS)
    \item 6: Mixture model (MIX)
    \item 7: BayesB and BayesA
    \item 8: genomic BLUP (GBLUP)
    \item 9: GWR is not applied. A normal prior distribution is assigned. The mean and variance are specified in Hyperpara.
    }
 See Details and References for these methods.
}
  \item{Referencevalues}{
  An Np length vector of reference values of parameters. These values are used to define the prior means and standard deviations of parameters at the first iteration in MPI. The prior means are defined as Referencevalues. The prior standard deviations are abs(Referencevalues) * 0.1. Referencevalues are also used to define the standard deviations of the proposed distributions for each parameter if SdforParameters is not specified. The standard deviations are also abs(Referencevalues) * 0.1. In addition, Referencevalues are used for checking Modelfunction. GenomeBasedModel evaluates Modelfunction using Input, Freevec, and Referencevalues, and stops if unusual values are output.
}
  \item{Modelfunction}{
  An R function or cpp function compiled by Rcpp that describes the model of interest. When PassMatrix is FALSE, this function should take three arguments, a column vector of Input, Freevec, and a vector of parameters with length Np, with this order and should output a vector of length D. When PassMatrix is TRUE, this function takes Input, Freevec, and a matrix of parameters for all genotypes as arguments, and outputs a matrix of D x Ny. This latter option is generally faster and allows parallel computation across genotypes within Modelfunction. See Experiment1 in the example for the details and examples.
}
  \item{PolygenicEffect}{
  An Np length logical vector to add a polygenic effect to genome-wide regression. Ignored when GBLUP is applied. The relationship matrix is created from Geno if K is not specified.
}
  \item{Transformation}{
  An Np length vector specifying parameter transformation methods in Modelfunction. If NULL, GenomeBasedModel considers that all parameters are not transformed. Elements of the vector are either of "nt" (not transformed), "log" (log transformation), or "logit" (logit transformation).}
  \item{Residualgroup}{
  A D length vector of integers starting from 1, indicating the groups of residual variances. The default is rep(1, D), which means that the residual variance is shared among the D observations.
  }
    \item{Hyperpara}{
  A nested list to specify the hyperparameters of regression methods in GWR for each parameter. The length of Hyperpara, length(Hyperpara), is Np and the length of the list for each parameter (e.g., length(Hyperpara[[1]]) for the first parameter) depends on the regression methods assumed. See Details for the details.}
  \item{Q}{
  An Nf x Ny matrix of covariates other than SNP genotypes in GWR (Nf is the number of covariates). If not specified, the intercept is automatically added.
  }
  \item{SdforParameters}{
  An Np length vector of the standard deviations of the proposed distributions for each parameter in MPI. These values influence the mixing of MCMC and may need to be tuned according to the acceptance rates shown on the console. The default is abs(Referencevalues) * 0.1.
  }
  \item{SdforVe}{
  An Nve length vector of the standard deviation of the proposed distribution for the residual variance in MPI. This value influences the mixing of MCMC and may need to be tuned according to the acceptance rate shown on the console. The default is abs(mean(Y)) * 0.05.
  }
  \item{K}{
  An (Nx + 1) x Nx matrix of the genetic relationship matrix used when GBLUP is assumed in GWR. THE FIRST LINE INCLUDES THE GENOTYPE IDs (only integers). If not specified although GBLUP is assumed, the relationship matrix is built from Geno using the A.mat function of the rrBLUP package.
}
  \item{Y_uncertainty}{
  A D x Ny matrix of uncertainties of observations. The order of genotypes should be same as Y. These values are added to the residual variance.
}
  \item{Nite_GWR}{
  An integer to specify the maximum number of iterations of GWR at one loop of MPI-GWR. See References for the convergence criterion.
}
  \item{Nloop}{
  An integer to specify the maximum number of MPI-GWR loops. The loop stops before Nloop if a convergence criterion is fulfilled (see References for convergence criterion).
}
  \item{Nite_MPI}{
  An integer to specify the number of iterations in MPI.
}
  \item{Nite_MPI_last}{
  An integer to specify the number of iterations in MPI, which is added after convergence. Samples obtained from these iterations are output as final results. Usually Nite_MPI_last > Nite_MPI.
}
  \item{Thin}{
  An integer to specify the number of sampling intervals in MPI. The number of samples is Nite_MPI/Thin during MPI-GWR loops and Nite_MPI_last/Thin in the final results.
}
  \item{Tol}{
  A scalar to specify the convergence criterion in GWR. See References for the details.
}
  \item{Nblock}{
  An integer to define the convergence criterion in MPI. Rhat is calculated from consecutive Nblock MPI sequences of likelihood. See References for Rhat and convergence criterion.
}
  \item{Rhatthr}{
  A scalar to define the convergence criterion in MPI. MPI-GWR loops stop if Rhat < Rhatthr for consecutive StopN loops. See References for the details.
}
  \item{StopN}{
  An integer to define the convergence criterion in MPI. MPI-GWR loops stop if Rhat < Rhatthr for consecutive StopN loops. See References for the details.
}
  \item{Lowerlimit}{
  An Np length vector to specify the lower limits of parameters. The default is -1e+20.
}
  \item{Upperlimit}{
  An Np length vector to specify the upper limits of parameters. The default is 1e+20.
}
  \item{Initialvalues}{
  A list object returned by GenomeBasedModel. If specified, the function restarts the MPI-GWR loops using the values in the list as the initial values. See Examples (Experiment 2) for an example.
}
  \item{PassMatrix}{
  A logical to determine the input and output data formats of Modelfunction. See Modelfunction.
}
}
\details{
  GenomeBasedModel repeats MPI and GWR loops until convergence or until the number of loops reaches Nloop-1. Then, the function conducts the final loop of MPI and GWR and outputs the results. Parameter inference is based on a variational Bayesian inference, but MPI is based on MCMC simulations. Eight regression methods are available for GWR. The details of model structures, parameter inference algorithms, and hyperparameters of these regression methods are provided in the vignette of VIGoR package (available from CRAN) or the online pdf of the package (https://github.com/Onogi/VIGoR). Note that wBSR, SSVS, and MIX are not supported by current VIGoR (ver. 1.1.0 later). \cr\cr
  Regression methods require different hyperparameters specified in Hyperpara. Below is the hyperparameters to be specified (default values) and simple examples are provided in Examples (Experiment 4).\cr
  \itemize{
    \item 1 (BL)    : Phi (1), Omega (1)
    \item 2 (EBL)   : Phi (0.1), Omega (0.1), Psi (1), Theta (0.1)
    \item 3 (wBSR)  : Nu (5), S2 (0.1), Kappa (0.01)
    \item 4 (BayesC): Nu (5), S2 (0.1), Kappa (0.01)
    \item 5 (SSVS)  : c (0.01), Nu (5), S2 (0.1), Kappa (0.01)
    \item 6 (MIX)   : c (0.01), Nu (5), S2 (0.1), Kappa (0.01)
    \item 7 (BayesB): Nu (5), S2 (0.1), Kappa (0.01)
    \item 8 (GBLUP) : Nu (-2), S2 (0)
    \item 9 (Normal prior) : Mu (0), Var (10000)
    }
  Default values are assigned to all hyperparameters and users can run the function without specifying hyperparameters. However, the default values do not ensure good performance of the model. Thus, users may need to tune the values, for example, using cross-validation. From my experience, BL, EBL, and BayesC work relatively well without tuning hyperparameters, and GBLUP is most insensitive to the choice of hyperparameters. When modifying the hyperparameters, note that the response variables (i.e., model parameters) are standardized in GWR.\cr
  Kappa indicates the probability that markers have non-zero effects. Thus, BayesB can be used as BayesA by setting Kappa to 1. Similarly, BayesC is Bayesian ridge regression when Kappa is 1.
}
\value{
  \item{$Para}{A length Np list of MCMC samples of model parameters in MPI. For each parameter, a Nite_MPI_last/Thin x Ny matrix is returned.}
  \item{$ResidualVar}{A length Nite_MPI_last/Thin vector of MCMC samples of residual variances in MPI.}
  \item{$Loglike}{A length Nite_MPI_last/Thin vector of MCMC samples of log likelihood in MPI}
  \item{$Genome}{A length Np list of estimated parameters in GWR. See the vignette of VIGoR package (available from CRAN) or the online pdf of the package (see references) for the details of the objects contained in the list. Generally, Beta denotes the posterior means of regression coefficients (SNP effects) and Alpha denotes the posterior means of coefficients for other covariates including the intercept. Sd.beta and Sd.alpha are the posterior standard deviations of these coefficients. Rho and Gamma denote the probability of each SNP to be included in the regression model. For GBLUP, U denotes the posterior means of genotypic values for each genotype and Sigma2 is the posterior mean of the variance component explained by U.}
  \item{$Final}{A length 3 list. This contains the means and variances of the model parameters used in the last MPI fitting, and the residual variance(s) at the last MPI iteration. These objects are used as initial values when given as Initialvalues. See Examples (Experiment 3) for an example.}
  \item{$Nite}{A Nloop x Np matrix storing the numbers of iterations conducted in GWR.}
}
\references{
  For model functions and convergence criterion,\cr
  Onogi A. Connecting mathematical models to genomes: Joint estimation of model parameters and genome-wide marker effects on these parameters. Bioinformatics 36, 3169-3176 (2020). \cr\cr
  For GWR, hyperparameters in GWR, and output objects by GWR,\cr
  Onogi A. Documents for VIGoR ver. 1.1.0, https://github.com/Onogi/VIGoR, (2021)\cr
}
\author{
  Akio Onogi
}
\examples{
#Basic settings#########################################################################
#Consider a sigmoid growth curve model which has parameters A and B.
#The function is
Sigmoid <- function (t, A, B) {B / (1 + exp(-A * t + 10))}
#at time t.

#The number of genotypes is 100
N <- 100


#Experiment 1###########################################################################
#Assume that every genotype has five observations measured at t=4, 8, 12, 16, and 20
Nsample <- 5

#Sampled time for each genotype
Time.1 <- matrix(seq(4, 20, 4), nrow = Nsample, ncol = N)

#Every genotype has 1000 SNP genotypes
P <- 1000
Geno <- matrix(sample (-1:1, P * N, replace = TRUE), P, N)

#Generate model parameters. For each parameter, 5 SNPs are randomly chosen as QTLs
Nqtl <- 5
QTL <- matrix(0, Nqtl, 2)
Effect <- matrix(0, Nqtl, 2)
for(parameter in 1:2){
  QTL[, parameter] <- sort(sample(1:1000, Nqtl, replace = FALSE))
  Effect[, parameter] <- rnorm(Nqtl, 0, 1)
}

A <- crossprod(Geno[QTL[, 1],], Effect[, 1, drop = FALSE])
B <- crossprod(Geno[QTL[, 2],], Effect[, 2, drop = FALSE])

#arrange the scales
A <- A/(max(A) - min(A)) * 1.5
A <- A - min(A) + 0.5
B <- B/(max(B) - min(B)) * 1.5
B <- B - min(B) + 0.5

#add noises. Heritability is 0.8
A.noise <- A + rnorm(N, 0, sqrt(var(A) * (1-0.8)/0.8))
B.noise <- B + rnorm(N, 0, sqrt(var(B) * (1-0.8)/0.8))

#create observations
Ve <- 0.01 #residual variance
Y.1 <- matrix(NA, Nsample, N)
for(n in 1:N){
  Y.1[1:Nsample, n] <- Sigmoid (Time.1[, n], A.noise[n], B.noise[n]) + rnorm (Nsample, 0, sqrt(Ve))
}

#plot(Time.1[, 1], Y.1[, 1], xlim = c(0, 20), ylim = c(0, 2.5), type = "o")
#for(n in 1:N){
#  points(Time.1[, n], Y.1[, n], type = "o")
#}

#create input objects for GenomeBasedModel
#Time.1 can be used as Input
Input <- Time.1

#Freevec can be used for arbitrary purposes, but not used here.
Freevec <- 0

#add IDs to Y and Geno
Y.1 <- rbind(1:N, Y.1)
Geno <- rbind(1:N, Geno)

#number of parameters (A and B)
Np <- 2

#Arbitrary numbers can be used as Missing (but NA should be avoided)
Missing <- 999999
#here Y.1 has no missing value

#Use EBL for parameter A and BayesB for B
Methodcode <- c(2, 7)

#Referencevalues is used to define the initial distributions
#and standard deviations of proposal distributions.
#this is also used for checking Modelfunction before analysis
#Enter not so bad values
Referencevalues <- c(1, 1)

#Model function
Model <- function(input, freevec, parameters){

  #input is a column vector of the argument Input

  a <- parameters[1]
  b <- parameters[2]

  output <- b / (1 + exp(-a * input + 10))

  #please ensure that output has length D.
  #and does not include NA and Inf.
  output
}

#run GenomeBasedModel
library(GenomeBasedModel)
\donttest{
#it will take a couple of minutes
Result.1 <- GenomeBasedModel(Input, Freevec, Y.1, Missing, Np,
Geno, Methodcode, Referencevalues, Model)

#compare with the true values
plot(A, apply(Result.1$Para[[1]], 2, mean));abline(0, 1)
plot(B, apply(Result.1$Para[[2]], 2, mean));abline(0, 1)

#SNP effects and QTL positions
plot(abs(Result.1$Genome[[1]]$Beta), main = "A", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 1], lty = 2)

plot(abs(Result.1$Genome[[2]]$Beta), main = "B", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 2], lty = 2)
#Rho of BayesB indicates the probability of model inclusion of SNPs
#and is often used for detecting QTLs
plot(abs(Result.1$Genome[[2]]$Rho), main = "B", xlab = "SNPs", ylab = "Rho (inclusion prob.)")
abline(v = QTL[, 2], lty = 2)
#Minor QTLs may be missed.

#Use the PassMatrix option
Model.PassMatrix <- function(input, freevec, parameters){

  #input is the argument Input
  #parameters is a Np x N matrix of parameter values
  n <- ncol(input)
  output <- NULL

  for(i in 1:n){
    a <- parameters[1, i]
    b <- parameters[2, i]
    output <- cbind(output, b / (1 + exp(-a * input[, i] + 10)))
  }

  output
}

Result.1.PassMatrix <- GenomeBasedModel(Input, Freevec, Y.1, Missing, Np,
Geno, Methodcode, Referencevalues, Model.PassMatrix, PassMatrix = TRUE)

#compare with the true values
plot(A, apply(Result.1.PassMatrix$Para[[1]], 2, mean));abline(0, 1)
plot(B, apply(Result.1.PassMatrix$Para[[2]], 2, mean));abline(0, 1)

#compare with Result.1
plot(apply(Result.1.PassMatrix$Para[[1]], 2, mean),
     apply(Result.1$Para[[1]], 2, mean));abline(0, 1)
plot(apply(Result.1.PassMatrix$Para[[2]], 2, mean),
     apply(Result.1$Para[[2]], 2, mean));abline(0, 1)
}

#Experiment 2###################################################################################
\dontrun{
#Assume that every genotype has three to five observations measured at arbitrary t (1 <= t <= 20)
Nsample <- sample(3:5, N, replace = TRUE)
Time.2 <- matrix(0, nrow = 5, ncol = N)
for(n in 1:N){
  Time.2[, n] <- c(sort(sample(1:20, Nsample[n], replace = FALSE)), rep(0, 5-Nsample[n]))
}

Y.2 <- matrix(0, 5, N)
for(n in 1:N){
  Y.2[, n] <- Sigmoid (Time.2[, n], A.noise[n], B.noise[n]) + rnorm (5, 0, sqrt(Ve))
}

#now Y.2 has missing values
Y.2[Time.2 == 0] <- Missing

Input <- Time.2
Y.2 <- rbind(1:N, Y.2)

#run GenomeBasedModel
Result.2 <- GenomeBasedModel(Input, Freevec, Y.2, Missing, Np,
Geno, Methodcode, Referencevalues, Model)

#Prolong iterations when the chain does not converge
Result.2 <- GenomeBasedModel(Input, Freevec, Y.2, Missing, Np,
Geno, Methodcode, Referencevalues, Model,
Initialvalues = Result.2)

#compare with the true values
plot(A, apply(Result.2$Para[[1]], 2, mean));abline(0, 1)
plot(B, apply(Result.2$Para[[2]], 2, mean));abline(0, 1)

#SNP effects and QTL positions
plot(abs(Result.2$Genome[[1]]$Beta), main = "A", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 1], lty = 2)

#plot(abs(Result.2$Genome[[2]]$Beta), main = "B", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 2], lty = 2)


#Experiment 3#####################################################################################
#try another format for Input and Y. Use the same data as experiment 2.
Time.3 <- matrix(1:20, nrow = 20, ncol = N)

Y.3 <- matrix(0, nrow = 20, ncol = N)
for(n in 1:N){
  Y.3[Time.2[, n], n] <- Y.2[-1, n][Y.2[-1, n] != Missing]
}
Y.3[Y.3 == 0] <- Missing

Input <- Time.3
Y.3 <- rbind(1:N, Y.3)

#run GenomeBasedModel
Result.3 <- GenomeBasedModel(Input, Freevec, Y.3, Missing, Np,
Geno, Methodcode, Referencevalues, Model)


#Experiment 4######################################################################################
#use other regression methods
Methodcode <- c(1, 4)#BL and BayesC

#change the hyperparameters
Hyperpara <- list (list(Phi = 1, Omega = 0.1), list(Nu=5, S2=0.1, Kappa=0.005))
#Bayesian ridge regression can be done by setting Kappa of BayesC to 1

Result.4 <- GenomeBasedModel(Input, Freevec, Y.3, Missing, Np,
Geno, Methodcode, Referencevalues, Model, Hyperpara = Hyperpara)

#compare with the true values
plot(A, apply(Result.4$Para[[1]], 2, mean));abline(0, 1)
plot(B, apply(Result.4$Para[[2]], 2, mean));abline(0, 1)

#SNP effects and QTL positions
plot(abs(Result.4$Genome[[1]]$Beta), main = "A", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 1], lty = 2)
#shrinkage of negative SNPs of BL is generally weaker than that of EBL.

plot(abs(Result.4$Genome[[2]]$Beta), main = "B", xlab = "SNPs", ylab = "Abs. effect size")
abline(v = QTL[, 2], lty = 2)
plot(abs(Result.4$Genome[[2]]$Rho), main = "B", xlab = "SNPs", ylab = "Inclusion probability")
abline(v = QTL[, 2], lty = 2)
#Inclusion probability (Rho) of BayesC and BayesB is often used as good signals of QTLs
}
}
