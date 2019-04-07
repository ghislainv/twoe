#===================================================================
# entry_ngrowth_gibbs.R
#
# entry_ngrowth_gibbs() samples from the posterior distribution of a
# Gaussian hierarchical linear regression model with observation errors
# in R using linked C++ code in Scythe.
#
# A hierarchical structure can be specififed.
# If so, the code uses Algorithm 2 of Chib & Carlin (1999) for efficient
# inference of (\beta | Y, sigma^2, Vb).
#
# Chib, S. & Carlin, B. P. (1999) On MCMC sampling in hierarchical
# longitudinal models, Statistics and Computing, 9, 17-26
#
#===================================================================
#
# Original code by Ghislain Vieilledent, June 2011
# CIRAD UR B&SEF
# ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
#
#===================================================================
# 
# This software is distributed under the terms of the GNU GENERAL
# PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
# file for more information.
#
# Copyright (C) 2011 Andrew D. Martin and Kevin M. Quinn
#
#===================================================================
#
# Revisions: 
# - G. Vieilledent, Sept 9 2011
#
#===================================================================


entry_ngrowth_gibbs <- function (fixed, random, group, diameter, data,
                                 burnin=1000, mcmc=10000, thin=10,
                                 verbose=1, seed=NA, a.sd=0.927,
                                 b.sd=0.0038, beta.start=NA,
                                 sigma2.start=NA, Vb.start=NA,
                                 mubeta=0, Vbeta=1.0E6, r, R,
                                 nu=0.001, delta=0.001, ...)  {
  #=====================  
  #= Fixed effects model
  if (is.null(random)) {
    
    #========
    # Basic checks
    #========
    check.mcmc.parameters.hmodels(burnin, mcmc, thin)
    check.verbose.hmodels(verbose)
    check.offset(list(...))

    #========
    # Seed
    #========
    seed <- form.seeds.hmodels(seed)
  
    #======== 
    # Form response and model matrices
    #========
    mf.fixed <- model.frame(formula=fixed,data=data)
    X <- model.matrix(attr(mf.fixed,"terms"),data=mf.fixed)
    Y <- model.response(mf.fixed)

    #======== 
    # Observation error process
    #========
    D <- check.diameter.hmodels(diameter,nrow(X))
    a.sd <- check.positive(a.sd)
    b.sd <- check.positive(b.sd)
    
    #======== 
    # Model parameters
    #========
    nobs <- nrow(X)
    np <- ncol(X)
    ngibbs <- mcmc+burnin
    nthin <- thin
    nburn <- burnin
    nsamp <- mcmc/thin

    #========
    # Form and check starting parameters
    #========
    beta.start <- form.beta.start.hmodels(fixed,data,beta.start,np,family="gaussian",defaults=NA)
    sigma2.start <- form.sigma2.start.hmodels(fixed,data,sigma2.start,family="gaussian")
      
    #========
    # Form priors
    #========
    mvn.prior <- form.mvn.prior.hmodels(mubeta,Vbeta,np)
    mubeta <- mvn.prior[[1]]
    Vbeta <- mvn.prior[[2]]
    check.ig.prior.hmodels(nu,delta)
    s1 <- nu
    s2 <- delta

    #========
    # Parameters to save
    #========
    beta_vect <- rep(c(beta.start),each=nsamp)
    V <- rep(sigma2.start,nsamp)
    Y_true <- rep(1,nobs) # !! must be >0 (cf. log(Y_true))
    log_Y_pred <- rep(0,nobs)
    Deviance <- rep(0,nsamp)

    #========
    # call C++ code to draw sample
    #========
    Sample <- .C("entry_ngrowth_gibbs_fixed",
                 #= Constants and data
                 ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn),## Number of iterations, burning and samples
                 nobs=as.integer(nobs), ## Constants
                 np=as.integer(np), ## Constants
                 Y_vect=as.double(c(Y)), ## Response variable
                 X_vect=as.double(c(X)), ## Covariates
                 D=as.double(c(D)), ## Diameter for observation errors
                 a_sd=as.double(a.sd), ## Intercept for observation errors
                 b_sd=as.double(b.sd), ## Slope for observation errors
                 #= Parameters to save
                 beta_vect.nonconst=as.double(beta_vect), ## Fixed parameters of the regression
                 V.nonconst=as.double(V), ## Variance of residuals
                 Y_true.nonconst=as.double(Y_true), ## Latent variable, conditional posterior mean
                 #= Defining priors
                 mubeta_vect=as.double(c(mubeta)), Vbeta_vect=as.double(c(Vbeta)),
                 s1_V=as.double(s1), s2_V=as.double(s2),
                 #= Diagnostic
                 Deviance.nonconst=as.double(Deviance),
                 log_Y_pred.nonconst=as.double(log_Y_pred), ## Predictive posterior mean
                 #= Seeds
                 seed=as.integer(seed), 
                 #= Verbose
                 verbose=as.integer(verbose),
                 PACKAGE="twoe")
 
    #= Matrix of MCMC samples
    Matrix <- matrix(NA,nrow=nsamp,ncol=np+2)
    names.fixed <- paste("beta.",colnames(X),sep="")
    names.variances <- "sigma2"
    colnames(Matrix) <- c(names.fixed,names.variances,"Deviance")

    #= Filling-in the matrix
    Matrix[,c(1:np)] <- matrix(Sample[[11]],ncol=np)
    Matrix[,ncol(Matrix)-1] <- Sample[[12]]
    Matrix[,ncol(Matrix)] <- Sample[[18]]

    #= Transform Sample list in an MCMC object
    MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)

    #= pD and DIC
    Y.true <- Sample[[13]]
    log.Y.pred <- Sample[[19]]
    Matrix.mean <- apply(Matrix,2,mean)
    beta.mean <- Matrix.mean[c(1:np)]
    V.mean <- Matrix.mean[ncol(Matrix)-1]
    logL.hat <- sum(dnorm(Y,Y.true,a.sd+b.sd*D,1))+sum(dnorm(log(Y.true),X%*%beta.mean,sd=sqrt(V.mean),1)) # log Likelihood
    D.hat <- -2*logL.hat
    D.bar <- Matrix.mean[ncol(Matrix)]
    pD <- as.numeric(D.bar-D.hat)
    DIC <- as.numeric(D.bar+pD)
  
    #= Output
    return (list(mcmc=MCMC,Y.true=Y.true,Y.pred=exp(log.Y.pred+0.5*V.mean),pD=pD,DIC=DIC))
  }

  #=====================
  #= Mixed effects model
  if (!is.null(random)) {
    
    #========
    # Basic checks
    #========
    check.group.hmodels(group, data)
    check.mcmc.parameters.hmodels(burnin, mcmc, thin)
    check.verbose.hmodels(verbose)
    check.offset(list(...))

    #========
    # Seed
    #========
    seed <- form.seeds.hmodels(seed)
  
    #======== 
    # Form response and model matrices
    #========
    mf.fixed <- model.frame(formula=fixed,data=data)
    X <- model.matrix(attr(mf.fixed,"terms"),data=mf.fixed)
    Y <- model.response(mf.fixed)
    mf.random <- model.frame(formula=random,data=data)
    W <- model.matrix(attr(mf.random,"terms"),data=mf.random)

    #======== 
    # Observation error process
    #========
    D <- check.diameter.hmodels(diameter,nrow(X))
    a.sd <- check.positive(a.sd)
    b.sd <- check.positive(b.sd)
    
    #======== 
    # Model parameters
    #========
    nobs <- nrow(X)
    IdentGroup <- as.numeric(as.factor(as.character(data[,names(data)==as.character(group)])))-1
    LevelsGroup <- sort(unique(IdentGroup+1))
    LevelsGroup.Name <- sort(unique(as.character(data[,names(data)==as.character(group)])))
    ngroup <- length(LevelsGroup)
    np <- ncol(X)
    nq <- ncol(W)
    ngibbs <- mcmc+burnin
    nthin <- thin
    nburn <- burnin
    nsamp <- mcmc/thin

    #========
    # Form and check starting parameters
    #========
    beta.start <- form.beta.start.hmodels(fixed,data,beta.start,np,family="gaussian",defaults=NA)
    sigma2.start <- form.sigma2.start.hmodels(fixed,data,sigma2.start,family="gaussian")
    Vb.start <- form.Vb.start.hmodels(Vb.start,nq)
      
    #========
    # Form priors
    #========
    mvn.prior <- form.mvn.prior.hmodels(mubeta,Vbeta,np)
    mubeta <- mvn.prior[[1]]
    Vbeta <- mvn.prior[[2]]
    wishart.prior <- form.wishart.prior.hmodels(r,R,nq)
    r <- wishart.prior[[1]]
    R <- wishart.prior[[2]]
    check.ig.prior.hmodels(nu,delta)
    s1 <- nu
    s2 <- delta

    #========
    # Parameters to save
    #========
    beta_vect <- rep(c(beta.start),each=nsamp)
    Vb_vect <- rep(c(Vb.start),each=nsamp)
    b_vect <- rep(0,nq*ngroup*nsamp)
    V <- rep(sigma2.start,nsamp)
    Y_true <- rep(1,nobs) # !! must be >0 (cf. log(Y_true))
    log_Y_pred <- rep(0,nobs)
    Deviance <- rep(0,nsamp)

    #========
    # call C++ code to draw sample
    #========
    Sample <- .C("entry_ngrowth_gibbs_mixed",
                 #= Constants and data
                 ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn),## Number of iterations, burning and samples
                 nobs=as.integer(nobs), ngroup=as.integer(ngroup), ## Constants
                 np=as.integer(np), nq=as.integer(nq), ## Constants
                 IdentGroup=as.integer(IdentGroup),
                 Y_vect=as.double(c(Y)), ## Response variable
                 X_vect=as.double(c(X)), ## Covariates
                 W_vect=as.double(c(W)), ## Covariates
                 D=as.double(c(D)), ## Diameter for observation errors
                 a_sd=as.double(a.sd), ## Intercept for observation errors
                 b_sd=as.double(b.sd), ## Slope for observation errors
                 #= Parameters to save
                 beta_vect.nonconst=as.double(beta_vect), ## Fixed parameters of the regression
                 b_vect.nonconst=as.double(b_vect), ## Random effects on intercept and slope
                 Vb_vect.nonconst=as.double(Vb_vect), ## Variance-covariance of random effects
                 V.nonconst=as.double(V), ## Variance of residuals
                 Y_true.nonconst=as.double(Y_true), ## Latent variable, conditional posterior mean
                 #= Defining priors
                 mubeta_vect=as.double(c(mubeta)), Vbeta_vect=as.double(c(Vbeta)),
                 r=as.double(r), R_vect=as.double(c(R)),
                 s1_V=as.double(s1), s2_V=as.double(s2),
                 #= Diagnostic
                 Deviance.nonconst=as.double(Deviance),
                 log_Y_pred.nonconst=as.double(log_Y_pred), ## Predictive posterior mean
                 #= Seeds
                 seed=as.integer(seed), 
                 #= Verbose
                 verbose=as.integer(verbose),
                 PACKAGE="twoe")
 
    #= Matrix of MCMC samples
    Matrix <- matrix(NA,nrow=nsamp,ncol=np+nq*ngroup+nq*nq+2)
    names.fixed <- paste("beta.",colnames(X),sep="")
    names.random <- paste("b.",rep(colnames(W),each=ngroup),".",rep(LevelsGroup.Name,nq),sep="")
    names.variances <- c(paste("VCV.",colnames(W),".",rep(colnames(W),each=nq),sep=""),"sigma2")
    colnames(Matrix) <- c(names.fixed,names.random,names.variances,"Deviance")

    #= Filling-in the matrix
    Matrix[,c(1:np)] <- matrix(Sample[[15]],ncol=np)
    Matrix[,c((np+1):(np+nq*ngroup))] <- matrix(Sample[[16]],ncol=nq*ngroup)
    Matrix[,c((np+nq*ngroup+1):(np+nq*ngroup+nq*nq))] <- matrix(Sample[[17]],ncol=nq*nq)
    Matrix[,ncol(Matrix)-1] <- Sample[[18]]
    Matrix[,ncol(Matrix)] <- Sample[[26]]

    #= Transform Sample list in an MCMC object
    MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)

    #= pD and DIC
    Y.true <- Sample[[19]]
    log.Y.pred <- Sample[[27]]
    Matrix.mean <- apply(Matrix,2,mean)
    beta.mean <- Matrix.mean[c(1:np)]
    bk.mean <- matrix(Matrix.mean[c((np+1):(np+nq*ngroup))],ncol=ngroup,nrow=nq,byrow=TRUE)
    V.mean <- Matrix.mean[ncol(Matrix)-1]
    logL.hat <- 0
    for (i in 1:nobs) {
      logL.hat <- logL.hat+
        dnorm(Y[i],Y.true[i],a.sd+b.sd*D[i],1)+
          dnorm(log(Y.true[i]),X[i,]%*%beta.mean+W[i,]%*%bk.mean[,c(IdentGroup[i]+1)],sd=sqrt(V.mean),1) # log Likelihood
    }
    D.hat <- -2*logL.hat
    D.bar <- Matrix.mean[ncol(Matrix)]
    pD <- as.numeric(D.bar-D.hat)
    DIC <- as.numeric(D.bar+pD)
  
    #= Output
    return (list(mcmc=MCMC,Y.true=Y.true,Y.pred=exp(log.Y.pred+0.5*V.mean),pD=pD,DIC=DIC))
  }
}

#===================================================================
# END
#===================================================================



