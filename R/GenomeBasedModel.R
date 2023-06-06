GenomeBasedModel <- function(Input, Freevec, Y, Missing, Np, Geno, Methodcode, Referencevalues, Modelfunction,
                             PolygenicEffect = FALSE, Transformation = NULL, Residualgroup = NULL, Hyperpara = NULL,
                             Q = NULL, SdforParameters = NULL, SdforVe = NULL, K = NULL, Y_uncertainty = NULL,
                             Nite_GWR = 1000, Nloop = 20, Nite_MPI = 500, Nite_MPI_last = 3000, Thin = 10, Tol = 1e-5,
                             Nblock = 5, Rhatthr = 1.01, StopN = 10,
                             Lowerlimit = rep(-1e+20, Np), Upperlimit = rep(1e+20, Np), Initialvalues = NULL, PassMatrix = FALSE){

  #check the arguments
  Ny <- ncol(Input)
  stopifnot(ncol(Input) == ncol(Y))
  LineID.Y <- Y[1, ]
  if(any(is.na(Input))) stop("Don't use NAs as missing values in Input")
  if(any(is.na(Y))) stop("Don't use NAs as missing values in Y")
  Nr <- nrow(Input)
  Ne <- nrow(Y) - 1
  stopifnot(Np == length(Referencevalues))
  stopifnot(Np == length(Methodcode))

  if(is.na(Missing)) stop ("Use an integer (e.g., 999999) for Missing values")

  if(!is.null(K)&!is.null(Geno)){
    stopifnot(identical(as.vector(Geno[1, ]), as.vector(K[1, ])))
    LineID.X <- K[1, ]
    K <- K[-1, ]
    Nx <- nrow(K)
  }else{
    if(is.null(K)&is.null(Geno)){
      stop("Specify either of Geno or K")
    }else{
      if(!is.null(K)){
        stopifnot(nrow(K[-1, ]) == ncol(K))
        LineID.X <- K[1, ]
        K <- K[-1, ]
        Nx <- nrow(K)
      }
      if(!is.null(Geno)){
        LineID.X <- Geno[1, ]
        Geno <- t(Geno[-1, ])
        Nx <- nrow(Geno)
      }
    }
  }

  Y_use <- match(LineID.Y, LineID.X)
  if(any(is.na(Y_use))) stop("Some IDs in Y does not exist in Geno/K")
  if(any(duplicated(LineID.Y))) stop("Some IDs in Y are duplicated")
  Y_use <- Y_use - 1
  Y <- Y[-1, ,drop=FALSE]

  X_use <- match(LineID.X, LineID.Y)
  if(any(is.na(X_use))) stop("Some IDs in Geno/K does not exist in Y")
  if(any(duplicated(LineID.X))) stop("Some IDs in Geno/K are duplicated")
  X_use <- X_use - 1

  #Polygenic effect
  PolygenicEffect <- as.logical(PolygenicEffect)
  if(length(PolygenicEffect) == 1){
    PolygenicEffect <- rep(PolygenicEffect, Np)
  } else {
    stopifnot(length(PolygenicEffect) == Np)
  }
  PolygenicEffect[Methodcode == 8] <- FALSE

  if(any(Methodcode == 8)|any(PolygenicEffect)){#Do GBLUP
    if(is.null(K)){
      Genotyperange <- range(Geno)
      if(Genotyperange[1] < (-1)|Genotyperange[2] > 1) stop("Genotype codes should be between -1 and 1 (required when calculating K)")
      cat("Create the genetic relationship matrix using A.mat\n")
      K <- A.mat(Geno, shrink = TRUE)
    }
    iK <- solve(K)
    iK.eigen <- eigen(iK)
  }else{
    iK <- 0
    iK.eigen <- list(values = 0, vectors = 0)
  }

  DoGWR <- which(Methodcode != 9)

  if(is.null(Geno)){
    Geno <- matrix(1, nrow = Nx, ncol = 1)
  }
  P <- ncol(Geno)

  #Check Modelfunction
  if(PassMatrix){
    Output <- Modelfunction(Input, Freevec, matrix(Referencevalues, nrow=Np, ncol=Ny))
  }else{
    Output <- matrix(0, nrow = Ne, ncol = Ny)
    for(line in 1:Ny){
      temp <- Modelfunction(Input[, line], Freevec, Referencevalues)
      if(length(temp) != Ne) stop("The length of Modelfunction output should be D (i.e., nrow(Y) - 1)")
      Output[, line] <- temp
    }
  }
  if(any(is.na(Output))) stop("Modelfunction outputs NA with Referencevalues")
  if(any(Output == Inf|Output == -Inf)) stop("Modelfunction outputs INf or -Inf with Referencevalues")

  #Transformation
  if(is.null(Transformation)){
    Transformation = rep(0, Np)
  }else{
    stopifnot(length(Transformation) == Np)
    if(any(!is.element(Transformation, c("nt","log","logit")))){
      stop("Elements in Transformation is either of nt, log, or logit")
    }
    Transformation [Transformation == "nt"] <- 0
    Transformation [Transformation == "log"] <- 1
    Transformation [Transformation == "logit"] <- 2
    Transformation <- as.numeric(Transformation)
  }

  #covariates for fixed effects
  if(is.null(Q)){
    Q <- rep(1, Ny)
    Q_x2 <- Ny
    Nf <- 1
  }else{
    stopifnot(ncol(Q) == ncol(Y))
    Q <- t(Q)
    Q_x2 <- colSums(Q^2)
    Nf <- ncol(Q)
  }

  #Homo/heteroscedasticity
  if(is.null(Residualgroup)){
    Nve <- 1
    Residualgroup <- rep(0, Ne)
  }else{
    if(length(Residualgroup) != Ne) stop("length of Residual group is different from the dimension of Y")
    if(length(unique(Residualgroup)) != max(Residualgroup)) stop("Use consecutive integers from 1 for Residualgroup")
    Nve <- max(Residualgroup)
    Residualgroup <- Residualgroup - 1
  }

  #hyperparameters for GWR
  if(is.null(Hyperpara)){
    #Default values
    Hyperpara <- as.list(numeric(Np))
    for(para in 1:Np){
      Hyperpara[[para]] <- switch(Methodcode[para],
                                  list(Phi = 1, Omega = 1, Psi = 0, Theta = 0, Nu = 0, S2 = 0, Kappa = 0, c = 0),#BL
                                  list(Phi = 0.1, Omega = 0.1, Psi = 1, Theta = 0.1, Nu = 0, S2 = 0, Kappa = 0, c = 0),#EBL
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu = 5, S2 = 0.1, Kappa = 0.01, c = 0),#wBSR
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu = 5, S2 = 0.1, Kappa = 0.01, c = 0),#BayesC
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu = 5, S2 = 0.1, Kappa = 0.01, c = 0.01),#SSVS
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu = 5, S2 = 0.1, Kappa = 0.01, c = 0.01),#MIX
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu = 5, S2 = 0.1, Kappa = 0.01, c = 0),#BayesB
                                  list(Phi = 0, Omega = 0, Psi = 0, Theta = 0, Nu= -2, S2 = 0, Kappa = 0, c = 0),#GBLUP
                                  list(Mu = 0, Var = 10000)#Normal prior
      )
    }
  }else{
    stopifnot(length(Hyperpara) == Np)
    for(para in 1:Np){
      switch(Methodcode[para],
             stopifnot(setequal(names(Hyperpara[[para]]), c("Omega", "Phi"))),#BL
             stopifnot(setequal(names(Hyperpara[[para]]), c("Omega", "Phi", "Theta", "Psi"))),#EBL
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2", "Kappa"))), #wBSR
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2", "Kappa"))), #BayesC
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2", "Kappa", "c"))), #SSVS
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2", "Kappa", "c"))), #MIX
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2", "Kappa"))), #BayesB
             stopifnot(setequal(names(Hyperpara[[para]]), c("Nu", "S2"))),#GBLUP
             stopifnot(setequal(names(Hyperpara[[para]]), c("Mu", "Var")))#Normal prior
      )
    }
    for(para in 1:Np){
      #add unused hyperparameters
      switch(Methodcode[para],
             Hyperpara[[para]] <- c(Hyperpara[[para]], Psi = 0, Theta = 0, Nu = 0, S2 = 0, Kappa = 0, c = 0),#BL
             Hyperpara[[para]] <- c(Hyperpara[[para]], Nu = 0, S2 = 0, Kappa = 0, c = 0),#EBL
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0, c = 0),#wBSR
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0, c = 0),#BayesC
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0),#SSVS
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0),#MIX
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0, c = 0),#BayesB
             Hyperpara[[para]] <- c(Hyperpara[[para]], Phi = 0, Omega = 0, Psi = 0, Theta = 0, Kappa = 0, c = 0)#GBLUP
      )
    }
  }

  #For scaling the model parameters variables in GWR
  Scalingvalues <- matrix(0, nrow = Np, ncol = 2)

  #Sd for proposal distributions in MPI
  if(is.null(SdforParameters)){
    SdforParameters <- 0.1 * abs(Referencevalues)
  }else{
    stopifnot(Np == length(SdforParameters))
  }

  if(is.null(SdforVe)){
    if(Nve == 1){
      SdforVe <- 0.05 * abs(mean(Y[Y!=Missing]))
    }else{
      temp <- Y
      temp[temp == Missing] <- NA
      SdforVe <- numeric(Nve)
      for(env in 1:Nve){
        SdforVe[env] <- 0.05 * abs(mean(temp[Residualgroup == env - 1, ], na.rm=TRUE))
      }
    }
  }else{
    if(Nve == 1){
      SdforVe <- SdforVe[1]
    }else{
      stopifnot(length(SdforVe) == length(unique(Residualgroup)))
    }
  }

  #check the number of iterations in MPI
  if(Nite_MPI > Nite_MPI_last){
    stop("Nite_MPI<Nite_MPI_last is recommended\n")
  }

  stopifnot(length(Lowerlimit)==Np)
  stopifnot(length(Upperlimit)==Np)

  #check Y_uncertainty
  if(is.null(Y_uncertainty)){
    Y_uncertainty <- matrix(0, nrow=Ne, ncol=Ny)
  }else{
    stopifnot(ncol(Y_uncertainty) == Ny)
    if(nrow(Y_uncertainty) != Ne) stop("The number of rows of Y_uncertainty should be nrow(Y) - 1")
  }

  #Objects used in MPI
  if(is.null(Initialvalues)){
    MeanOfPrior <- matrix(Referencevalues, nrow = Np, ncol = Ny)
    VarOfPrior <- 0.1 * abs(Referencevalues)
    Ve <- rep(0.1, Nve)
  }else{
    MeanOfPrior <- Initialvalues$Final$GWRfitting
    VarOfPrior <- Initialvalues$Final$GWRresidual
    Ve <- Initialvalues$Final$MPIresidual
  }
  if(any(Methodcode == 9)){
    Which <- which(Methodcode == 9)
    for(para in Which){
      MeanOfPrior[para, ] <- Hyperpara[[para]]$Mu
      VarOfPrior[para] <- Hyperpara[[para]]$Var
    }
  }

  Parameters <- matrix(0, nrow = Np, ncol = Ny)
  Ns <- Nite_MPI/Thin
  Ns_last <- Nite_MPI_last/Thin
  Loglike <- numeric(max(Ns, Ns_last))
  SampledPara <- matrix(0, nrow = Np * length(Loglike), ncol = Ny)
  SampledVe <- matrix(0, nrow = length(Loglike), ncol = Nve)
  Burnin <- 0#Ns and Ns_last should be changed when Burnin is not 0.
  Order.Ny <- 1:Ny
  Order.Ny <- Order.Ny-1#indices start from 0 in C
  NeNl <- numeric(Nve)
  v <- rowSums(Y!=Missing)
  for(env in 1:Nve){
    NeNl[env]=sum(v[Residualgroup == env-1])#indices in Residualgroup start from 0
  }

  #For assessing convergence
  Psij <- numeric(Nblock)
  Sj2 <- numeric(Nblock)
  Rhat <- rep(1e+8, StopN)
  Converge <- FALSE

  #objects used in GWR
  CondResidual <- 1
  Priortype <- numeric(Np)
  for(para in DoGWR){
    switch(Methodcode[para],
           Priortype[para] <- 1,
           Priortype[para] <- 1,
           Priortype[para] <- 2,
           Priortype[para] <- 2,
           Priortype[para] <- 2,
           Priortype[para] <- 2,
           Priortype[para] <- 2,
           Priortype[para] <- 3)
  }
  Y_variance <- matrix(0, nrow = Np, ncol = Ny)
  Y_expErrors <- numeric(Ny)
  X_x2 <- colSums(Geno^2)
  Result.GWR <- as.list(numeric(Np))
  X_expEffect <- as.list(numeric(Np))
  X_varEffect <- as.list(numeric(Np))
  X_exp2Effect <- as.list(numeric(Np))
  X_expGamma <- as.list(numeric(Np))
  X_exp2Gamma <- as.list(numeric(Np))
  X_expTau2 <- as.list(numeric(Np))
  X_expInTau2 <- as.list(numeric(Np))
  X_expEta2 <- as.list(numeric(Np))
  X_expSigma2 <- as.list(numeric(Np))
  X_s2 <- as.list(numeric(Np))
  Q_expEffect <- as.list(numeric(Np))
  Q_varEffect <- as.list(numeric(Np))
  Q_exp2Effect <- as.list(numeric(Np))
  Tau0 <- as.list(numeric(Np))
  expDelta2 <- as.list(numeric(Np))
  U_expBv <- as.list(numeric(Np))
  U_varBv <- as.list(numeric(Np))
  U_expSigma2 <- as.list(numeric(Np))
  U_varSigma2 <- as.list(numeric(Np))

  for(para in DoGWR){
    switch(Priortype[para],
           X_expEffect[[para]] <- numeric(P),
           X_expEffect[[para]] <- numeric(P),
           X_expEffect[[para]] <- 0
    )
    switch(Priortype[para],
           X_varEffect[[para]] <- numeric(P),
           X_varEffect[[para]] <- numeric(P),
           X_varEffect[[para]] <- 0
    )
    switch(Priortype[para],
           X_exp2Effect[[para]] <- numeric(P),
           X_exp2Effect[[para]] <- numeric(P),
           X_exp2Effect[[para]] <- 0
    )
    switch(Priortype[para],
           X_expGamma[[para]] <- 0,
           X_expGamma[[para]] <- rep(1, P),
           X_expGamma[[para]] <- 0
    )
    switch(Priortype[para],
           X_exp2Gamma[[para]] <- 0,
           X_exp2Gamma[[para]] <- rep(1, P),
           X_exp2Gamma[[para]] <- 0
    )
    switch(Priortype[para],
           X_expTau2[[para]] <- rep(P, P),
           X_expTau2[[para]] <- 0,
           X_expTau2[[para]] <- 0
    )
    switch(Priortype[para],
           X_expInTau2[[para]] <- rep(1/P, P),
           X_expInTau2[[para]] <- 0,
           X_expInTau2[[para]] <- 0
    )
    switch(Priortype[para],
           X_expEta2[[para]] <- rep(1, P),
           X_expEta2[[para]] <- 0,
           X_expEta2[[para]] <- 0
    )
    switch(Priortype[para],
           X_expSigma2[[para]] <- 0,
           X_expSigma2[[para]] <- rep(1, P),
           X_expSigma2[[para]] <- 0
    )
    switch(Priortype[para],
           X_s2[[para]] <- 0,
           X_s2[[para]] <- rep(1, P),
           X_s2[[para]] <- 0
    )

    Q_expEffect[[para]] <- numeric(Nf)
    Q_varEffect[[para]] <- numeric(Nf)
    Q_exp2Effect[[para]] <- numeric(Nf)
    Tau0[[para]] <- c(100, 1)
    expDelta2[[para]] <- 1

    if(PolygenicEffect[para]){
      U_expBv[[para]] <- rep(0, Nx)
      U_varBv[[para]] <- as.vector(diag(Nx))
    }else{
      switch(Priortype[para],
             U_expBv[[para]] <- 0,
             U_expBv[[para]] <- 0,
             U_expBv[[para]] <- rep(0, Nx)
      )
      switch(Priortype[para],
             U_varBv[[para]] <- 0,
             U_varBv[[para]] <- 0,
             U_varBv[[para]] <- as.vector(diag(Nx))
      )
    }

    U_expSigma2[[para]] <- 1
    U_varSigma2[[para]] <- 1
  }

  Order.P <- 1:P
  Order.P <- Order.P-1#indices start from 0 in C

  Nite_GWR.summary <- matrix(0, nrow=Nloop, ncol=Np)

  #Start GWR-MPI loops
  loop <- 0
  repeat{
    loop <- loop + 1
    if(loop<Nloop) {
      cat("Loop",  loop, "\n")
      Ns_loop <- Ns
      Nite_MPI_loop <- Nite_MPI
    }else{
      cat("Last loop\n")
      Ns_loop <- Ns_last
      Nite_MPI_loop <- Nite_MPI_last
    }

    #MPI
    cat("MPI\n")
    AccParameters <- numeric(Np)
    AccVe <- numeric(Nve)

    nRandom_para <- matrix(rnorm((Nite_MPI_loop + 1) * Np * Ny), nrow = Nite_MPI_loop + 1)
    nRandom_ve <- matrix(rnorm(Nite_MPI_loop * Nve), nrow = Nite_MPI_loop)
    uRandom_para <- matrix(runif(Nite_MPI_loop * Np * Ny), nrow = Nite_MPI_loop)
    uRandom_ve <- matrix(runif(Nite_MPI_loop * Nve), nrow = Nite_MPI_loop)

    if(PassMatrix){
      Result.MPI <- ModelParameterInference_PassMatrix(Input, Freevec, Y, Missing, Y_uncertainty,
                                                       Ny, Ne, Np, Nr,
                                                       SampledPara, SampledVe,
                                                       MeanOfPrior, VarOfPrior, Ve, Nve, Residualgroup,
                                                       SdforParameters, SdforVe, Loglike,
                                                       Burnin, Thin, Nite_MPI_loop, Order.Ny,
                                                       Lowerlimit, Upperlimit, AccParameters, AccVe, Modelfunction,
                                                       NeNl, Transformation,
                                                       nRandom_para, nRandom_ve, uRandom_para, uRandom_ve)
    }else{
      Result.MPI <- ModelParameterInference(Input, Freevec, Y, Missing, Y_uncertainty,
                                            Ny, Ne, Np, Nr,
                                            SampledPara, SampledVe,
                                            MeanOfPrior, VarOfPrior, Ve, Nve, Residualgroup,
                                            SdforParameters, SdforVe, Loglike,
                                            Burnin, Thin, Nite_MPI_loop, Order.Ny,
                                            Lowerlimit, Upperlimit, AccParameters, AccVe, Modelfunction,
                                            NeNl, Transformation,
                                            nRandom_para, nRandom_ve, uRandom_para, uRandom_ve)
    }

    if(loop < Nloop){
      Psij <- c(Psij, mean(Loglike[1:Ns]))
      Psij <- Psij[-1]
      Sj2 <- c(Sj2, sum((Loglike[1:Ns]-Psij[Nblock])^2)/(Ns-1))
      Sj2 <- Sj2[-1]
      if(loop >= Nblock){
        B <- sum((Psij - mean(Psij))^2)/(Nblock-1) * Ns
        W <- mean(Sj2)
        Rhat <- c(Rhat, sqrt(((Ns-1)/Ns * W + B/Ns)/W))
        Rhat <- Rhat[-1]
        cat("Rhat:", Rhat[StopN], "\n")
      }
    }
    cat("Mean loglike:", mean(Loglike[1:Ns_loop]), "\n")
    temp <- SampledPara[1:(Np * Ns_loop), ,drop=FALSE]
    cat("Residual var:",  Ve,  "\n")
    cat("Acceptance rate (parameters):",  AccParameters/(Ny * Nite_MPI_loop), "\n")
    cat("Acceptance rate (residual var):",  AccVe/(Nite_MPI_loop), "\n")
    for(para in 1:Np){
      Parameters[para, ] <- apply(temp[((para-1) * Ns_loop + 1):(para * Ns_loop), ,drop=FALSE], 2, mean)
      Y_variance[para, ] <- apply(temp[((para-1) * Ns_loop + 1):(para * Ns_loop), ,drop=FALSE], 2, var)
    }

    #GWR
    for(para in DoGWR){
      cat("GWR para", para)
      Scalingvalues[para, 1] <- mean(Parameters[para, ])
      Scalingvalues[para, 2] <- sd(Parameters[para, ])
      Y_stobs <- (Parameters[para, ]-Scalingvalues[para, 1])/Scalingvalues[para, 2]
      Y_variance[para, ] <- Y_variance[para, ]/Scalingvalues[para, 2]^2
      Niterations <- c(Nite_GWR, 0)
      Result.GWR[[para]] <- .C("GenomeWideRegression",
                             as.integer(Priortype[para]),
                             as.integer(Methodcode[para]),
                             as.integer(CondResidual),
                             as.integer(P),
                             as.integer(Nf),
                             as.integer(Ny),
                             as.integer(Nx),
                             as.integer(Niterations),#8
                             as.double(Y_stobs),
                             as.integer(Y_use),
                             as.double(Y_variance[para, ]),
                             as.double(Y_expErrors),#12. Initialized in the C function.
                             as.double(Geno),
                             as.double(X_x2),
                             as.double(X_expEffect[[para]]),#15
                             as.double(X_varEffect[[para]]),#16
                             as.double(X_exp2Effect[[para]]),#17
                             as.double(X_expGamma[[para]]),#18
                             as.double(X_exp2Gamma[[para]]),#19
                             as.double(X_expTau2[[para]]),#20
                             as.double(X_expInTau2[[para]]),#21
                             as.double(X_expEta2[[para]]),#22
                             as.double(X_expSigma2[[para]]),#23
                             as.double(X_s2[[para]]),#24
                             as.double(Q),
                             as.double(Q_expEffect[[para]]),#26
                             as.double(Q_varEffect[[para]]),#27
                             as.double(Q_exp2Effect[[para]]),#28
                             as.double(Q_x2),
                             as.double(Hyperpara[[para]]$Phi),
                             as.double(Hyperpara[[para]]$Omega),
                             as.double(Hyperpara[[para]]$Psi),
                             as.double(Hyperpara[[para]]$Theta),
                             as.double(Hyperpara[[para]]$Nu),
                             as.double(Hyperpara[[para]]$S2),
                             as.double(Hyperpara[[para]]$Kappa),
                             as.double(Hyperpara[[para]]$c),
                             as.double(expDelta2[[para]]),#38
                             as.double(Tau0[[para]]),#39
                             as.integer(Order.P),
                             as.double(iK.eigen$values),
                             as.double(iK.eigen$vectors),
                             as.double(t(iK.eigen$vectors)),
                             as.double(U_expBv[[para]]),#44
                             as.double(U_varBv[[para]]),#45
                             as.integer(X_use),
                             as.double(U_expSigma2[[para]]),#47
                             as.double(U_varSigma2[[para]]),#48
                             as.double(Tol),
                             as.integer(PolygenicEffect[para]),
                             PACKAGE = "GenomeBasedModel"
      )

      #update
      X_expEffect[[para]] <- Result.GWR[[para]][[15]]
      X_varEffect[[para]] <- Result.GWR[[para]][[16]]
      X_exp2Effect[[para]] <- Result.GWR[[para]][[17]]
      X_expGamma[[para]] <- Result.GWR[[para]][[18]]
      X_exp2Gamma[[para]] <- Result.GWR[[para]][[19]]
      X_expTau2[[para]] <- Result.GWR[[para]][[20]]
      X_expInTau2[[para]] <- Result.GWR[[para]][[21]]
      X_expEta2[[para]] <- Result.GWR[[para]][[22]]
      X_expSigma2[[para]] <- Result.GWR[[para]][[23]]
      X_s2[[para]] <- Result.GWR[[para]][[24]]
      Q_expEffect[[para]] <- Result.GWR[[para]][[26]]
      Q_varEffect[[para]] <- Result.GWR[[para]][[27]]
      Q_exp2Effect[[para]] <- Result.GWR[[para]][[28]]
      expDelta2[[para]] <- Result.GWR[[para]][[38]]
      Tau0[[para]] <- Result.GWR[[para]][[39]]
      U_expBv[[para]] <- Result.GWR[[para]][[44]]
      U_varBv[[para]] <- Result.GWR[[para]][[45]]
      U_expSigma2[[para]] <- Result.GWR[[para]][[47]]
      U_varSigma2[[para]] <- Result.GWR[[para]][[48]]

      cat(" residualvar (scaled):",  round(1/Tau0[[para]][1], 5), " Nite:", Result.GWR[[para]][[8]][2], "\n")
      Nite_GWR.summary[loop, para] <- Result.GWR[[para]][[8]][2]

      MeanOfPrior[para, ] <- (Y_stobs-Result.GWR[[para]][[12]]) * Scalingvalues[para, 2] + Scalingvalues[para, 1]
      VarOfPrior[para] <- 1/Tau0[[para]][1] * Scalingvalues[para, 2]^2
    }#para
    cat("\n")
    if(loop == Nloop) break
    if(!any(Rhat >= Rhatthr)) {
      loop <- Nloop-1
      Converge <- TRUE
    }
  }#loop

  List.markers <- as.list(numeric(Np))#scales are not backed
  for(para in DoGWR){
    switch(Methodcode[para],
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#BL
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Tau2 = Result.GWR[[para]][[20]],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      Lambda2 = Result.GWR[[para]][[38]],
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#EBL
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Tau2 = Result.GWR[[para]][[20]],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      Delta2 = Result.GWR[[para]][[38]],
                                      Eta2 = Result.GWR[[para]][[22]],
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]] * Result.GWR[[para]][[18]],#wBSR
                                      Sd.beta = sqrt(Result.GWR[[para]][[18]] * (1-Result.GWR[[para]][[18]]) * Result.GWR[[para]][[17]] +
                                                     Result.GWR[[para]][[18]]^2 * Result.GWR[[para]][[16]]),
                                      Gamma = Result.GWR[[para]][[18]],
                                      Sigma2 = Result.GWR[[para]][[23]],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ReidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#BayesC
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Rho = Result.GWR[[para]][[18]],
                                      Sigma2 = Result.GWR[[para]][[23]][1],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#SSVS
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Rho = Result.GWR[[para]][[18]],
                                      Sigma2 = Result.GWR[[para]][[23]][1],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#MIX
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Rho = Result.GWR[[para]][[18]],
                                      Sigma2 = Result.GWR[[para]][[23]][1:2],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(Beta = Result.GWR[[para]][[15]],#BayesB
                                      Sd.beta = sqrt(Result.GWR[[para]][[16]]),
                                      Rho = Result.GWR[[para]][[18]],
                                      Sigma2 = Result.GWR[[para]][[23]],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1]),
           List.markers[[para]] <- list(U = Result.GWR[[para]][[44]],#GBLUP
                                      Sd.u = sqrt(diag(matrix(Result.GWR[[para]][[45]],nrow = Nx))),
                                      Sigma2.u = Result.GWR[[para]][[47]],
                                      Alpha = Result.GWR[[para]][[26]],
                                      Sd.alpha = sqrt(Result.GWR[[para]][[27]]),
                                      ResidualVar = 1/Result.GWR[[para]][[39]][1])
    )
    if(PolygenicEffect[para]){
      List.markers[[para]] <- c(List.markers[[para]],
                                list(U = Result.GWR[[para]][[44]],
                                     Sd.u = sqrt(diag(matrix(Result.GWR[[para]][[45]],nrow = Nx))),
                                     Sigma2.u = Result.GWR[[para]][[47]])
                                )
    }
  }
  SampledPara <- matrix(SampledPara, ncol = Ny)
  List.parameters <- as.list(numeric(Np))
  for(para in 1:Np){
    List.parameters[[para]] <- SampledPara[((para-1) * Ns_last + 1):(para * Ns_last),]
  }

  Result <- list(Para = List.parameters,
                 ResidualVar = SampledVe,
                 Loglike = Loglike,
                 Genome = List.markers,
                 Final = list(GWRfitting = MeanOfPrior,GWRresidual = VarOfPrior,MPIresidual = Ve),
                 Nite = Nite_GWR.summary
  )
  rm(List.markers,Result.GWR,Result.MPI)
  rm(nRandom_para, nRandom_ve, uRandom_para, uRandom_ve)

  cat("=======================================\n")
  if(Converge){
    cat("Rhat in consective", StopN, "loops suggest the convergence of likelihood\n")
    cat("Rhat\n")
    cat(Rhat, "\n")
  }else{
    cat("Rhat in consective", StopN, "loops does not suggest the convergence of likelihood\n")
    cat("Rhat\n")
    cat(Rhat, "\n")
    cat("The process terminated because the number of loops reached", Nloop, "\n")
    cat("Recommend continueing calculation by giving the result list to GenomeBasedModel\n")
    cat("as the argument Initialvalues\n")
  }
  cat("=======================================\n")

  return(Result)
}
