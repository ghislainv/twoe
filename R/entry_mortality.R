entry_mortality <- function (data_mortality,burnin=5000,mcmc=5000,thin=5,th.sign=5) {
   
  # Some necessary objects
  DataMortality <- data_mortality
  DataMortality$IdentSpecies <- as.numeric(as.factor(as.character(DataMortality$Sp)))
  Levels.Species <- sort(unique(DataMortality$IdentSpecies))
  nspecies <- length(Levels.Species)
  Levels.Species.Name <- sort(unique(as.character(DataMortality$Sp)))

  # Transformed variables
  D <- DataMortality$D_t-20
  C <- DataMortality$C_t-20

  ###############
  # Gibbs sampler
  
  # Y, Status dead (1) or alive (0)
  # X1, (DataMortality$D_t-20) # cm
  # X2, (DataMortality$C_t-20) # m2.ha-1
   
  Start.C.Process <- Sys.time()
  model <-  entry_mortality_gibbs(fixed=status_tp1~D+C,
                                  random=~D+C,
                                  group="Sp",
                                  interval=DataMortality$interval_tp1,
                                  data=DataMortality,
                                  burnin=burnin, mcmc=mcmc, thin=thin,
                                  verbose=1, seed=NA,
                                  beta.start=c(0,0,0), sigma2.start=1, Vb.start=diag(1,3), mubeta=c(0,0,0),
                                  Vbeta=diag(1.0E6,3), r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001,
                                  FixOD=1)
  End.C.Process <- Sys.time()

  ########################
  # Computation time
  sink(file="time_mortality.txt")
  Time.C.Process <- difftime(End.C.Process,Start.C.Process,units="mins")
  cat("N# iterations: ",burnin+mcmc,"\n","Duration of the C++ program (in mins): ",Time.C.Process,sep="")
  sink()
  
  ##################################################################################################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #                                     Analysing MCMC ouputs                                      #
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ################################################################################################

  cat("\nAnalysing MCMC outputs and plotting diagnostic graphs...\n")
  flush.console()

  ########################
  # Change the parameter names for the mcmc object
  MCMC <- as.data.frame(model$mcmc)
  names.fixed <- paste("beta",c(0:2),sep="")
  names.random <- paste("b",rep(c(0:2),each=nspecies),".",rep(Levels.Species,3),sep="")
  names.variances <- c(paste("VCV.",c(0:2),rep(c(0:2),each=3),sep=""),"sigma2")
  colnames(MCMC) <- c(names.fixed,names.random,names.variances,"Deviance")
  mcmc.obj <- mcmc(MCMC,start=burnin+1,end=mcmc+burnin,thin=thin)
  
  ########################
  # Graphics and summaries
  pdf(file="posteriors_mortality.pdf")
  plot(mcmc.obj)
  dev.off()

  ########################
  # Export MCMCs
  write.table(MCMC,file="MCMC_mortality.txt",sep="\t",row.names=FALSE)
  # MCMC <- read.table(file="MCMC_mortality.txt",sep="\t",header=TRUE)
 
  ########################################################################################################################
  ########################################################################################################################
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##
  ##                                             Graphics
  ##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ########################################################################################################################
  ########################################################################################################################

  # Maximum likelihood approach
  Binom.Likelihood <- function (mu,Y.dead,Y.alive) {
    logvrais <- sum(log(1-(1-mu)^Y.dead))+sum(log((1-mu)^Y.alive))
    if (is.na(logvrais)) {logvrais <- -1.0E6}
    else {
      if (logvrais==Inf) {logvrais <- 1.0E6}
      if (logvrais==-Inf) {logvrais <- -1.0E6}
    } 
    return(logvrais)
  }

  # Correction factor
  # See Diggle P, Heagerty P, Liang K, Zeger S (2004). Analysis of Longitudinal Data.
  # and Hadfield's course notes, page 46.
  # http://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf
  # c2 <- (16 * sqrt(3)/(15 * pi))^2
  # CF <- 1/sqrt(1 + c2 * MCMC$V)
  # with MCMC$V=1 => CF <- 1/sqrt(1+(16*sqrt(3)/(15*pi))^2)
  CF_mortality <- 1/sqrt(1+(16*sqrt(3)/(15*pi))^2)
  
  # D
  pdf("theta_dbh.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))

  for (k in 1:nspecies) {
    
    # C fixed to the mean of the population
    C.fixed <- round(mean(DataMortality$C_t[DataMortality$IdentSpecies==Levels.Species[k]]))

    # Row data
    D.sp <- c(DataMortality$D_t[DataMortality$IdentSpecies==Levels.Species[k]])
    i.sp <- c(DataMortality$interval_tp1[DataMortality$IdentSpecies==Levels.Species[k]])
    s.sp <- c(DataMortality$status_tp1[DataMortality$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(D.sp); xmax.Data <- max(D.sp)

    # Bins for MLL mean mortality rate by DBH class
    bins <- c(min(10,D.sp-1),20,40,60,100,max(200,D.sp+1))
    midpoints <- c(min(10,D.sp-1)+(20-min(10,D.sp-1))/2,30,50,80,100+(max(200,D.sp+1)-100)/2)
    dbhsplit <- cut(D.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    levels.bin <- levels(dbhsplit)
    MLL.theta <- rep(0,length(midpoints))
    for (i in 1:length(midpoints)) {
      Y.dead <- i.sp[dbhsplit==levels.bin[i] & s.sp==1]
      Y.alive <- i.sp[dbhsplit==levels.bin[i] & s.sp==0]
      MinMLL <- optimize(f=Binom.Likelihood,lower=1.0E-6,upper=1,Y.dead=Y.dead,Y.alive=Y.alive,maximum=TRUE)
      if (length(Y.dead)==0&length(Y.alive)==0) {MLL.theta[i] <- NA}
      if (length(Y.dead)==0&length(Y.alive)!=0) {MLL.theta[i] <- 0}
      if (length(Y.dead)!=0&length(Y.alive)!=0) {MLL.theta[i] <- MLL.theta[i] <- MinMLL$maximum}
    }
    ymin.MLL <- min(MLL.theta,na.rm=TRUE); ymax.MLL <- max(MLL.theta,na.rm=TRUE);
    index1 <- which(!is.na(MLL.theta))[1]
    index2 <- which(!is.na(MLL.theta))[length(which(!is.na(MLL.theta)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    D.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)

    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")])
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    theta.hat <- vector()
    for (n in 1:100) {
      theta.hat[n] <-  inv_logit((b0sp+b1sp*(D.seq[n]-20)+b2sp*(C.fixed-20))*CF_mortality)
    }
    ymin.Pred <- min(theta.hat); ymax.Pred <- max(theta.hat)

    # Prediction mean species
    b0.msp <- mean(MCMC$beta0)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    theta.hat.msp <- vector()
    for (n in 1:100) {
      theta.hat.msp[n] <-  inv_logit((b0.msp+b1.msp*(D.seq[n]-20)+b2.msp*(C.fixed-20))*CF_mortality)
    }
    ymin.Pred.msp <- min(theta.hat.msp); ymax.Pred.msp <- max(theta.hat.msp)
    ymin.Graph <- min(ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.MLL,ymax.Pred,ymax.Pred.msp)

    plot(midpoints,MLL.theta,
         pch=18,
         cex=0.8,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="DBH (cm)",
         ylab="Annual mortality rate")
    lines(D.seq,theta.hat.msp,col="grey")
    lines(D.seq,theta.hat)
    text(x=xmin.Graph,
         y=ymin.Graph+(ymax.Graph-ymin.Graph)*0.90,
         labels=paste("C mean=",C.fixed," m2.ha-1",sep=""),cex=0.8,pos=4)   
  }
  dev.off()

  # C
  pdf("theta_comp.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))
  
  for (k in 1:nspecies) {
    
    # D fixed to the mean of the population
    D.fixed <- round(mean(DataMortality$D_t[DataMortality$IdentSpecies==Levels.Species[k]]))

    # Row data
    C.sp <- c(DataMortality$C_t[DataMortality$IdentSpecies==Levels.Species[k]])
    i.sp <- c(DataMortality$interval_tp1[DataMortality$IdentSpecies==Levels.Species[k]])
    s.sp <- c(DataMortality$status_tp1[DataMortality$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(C.sp); xmax.Data <- max(C.sp)
    
    # Bins for MLL mean mortality rate by Comp class
    bins <- c(0,10,20,30,50,70,max(100,C.sp+1))
    midpoints <- c(5,15,25,45,60,70+(max(100,C.sp+1)-70)/2)
    Csplit <-  cut(C.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    levels.bin <- levels(Csplit)
    MLL.theta <- rep(0,length(midpoints))
    for (i in 1:length(midpoints)) {
      Y.dead <- i.sp[Csplit==levels.bin[i] & s.sp==1]
      Y.alive <- i.sp[Csplit==levels.bin[i] & s.sp==0]
      MinMLL <- optimize(f=Binom.Likelihood,lower=1.0E-6,upper=1,Y.dead=Y.dead,Y.alive=Y.alive,maximum=TRUE)
      if (length(Y.dead)==0&length(Y.alive)==0) {MLL.theta[i] <- NA}
      if (length(Y.dead)==0&length(Y.alive)!=0) {MLL.theta[i] <- 0}
      if (length(Y.dead)!=0&length(Y.alive)!=0) {MLL.theta[i] <- MLL.theta[i] <- MinMLL$maximum}
    }
    ymin.MLL <- min(MLL.theta,na.rm=TRUE); ymax.MLL <- max(MLL.theta,na.rm=TRUE);
    index1 <- which(!is.na(MLL.theta))[1]
    index2 <- which(!is.na(MLL.theta))[length(which(!is.na(MLL.theta)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    C.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)
 
    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")])
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    theta.hat <- vector()
    for (n in 1:100) {
      theta.hat[n] <-  inv_logit((b0sp+b1sp*(D.fixed-20)+b2sp*(C.seq[n]-20))*CF_mortality)
    }
    ymin.Pred <- min(theta.hat); ymax.Pred <- max(theta.hat)

    # Prediction mean species
    b0.msp <- mean(MCMC$beta0)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    theta.hat.msp <- vector()
    for (n in 1:100) {
      theta.hat.msp[n] <-  inv_logit((b0.msp+b1.msp*(D.fixed-20)+b2.msp*(C.seq[n]-20))*CF_mortality)
    }
    ymin.Pred.msp <- min(theta.hat.msp); ymax.Pred.msp <- max(theta.hat.msp)
    ymin.Graph <- min(ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.MLL,ymax.Pred,ymax.Pred.msp)

    plot(midpoints,MLL.theta,
         pch=18,
         cex=0.8,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="Competition (m2.ha-1)",
         ylab="Annual mortality rate")
    lines(C.seq,theta.hat.msp,col="grey")
    lines(C.seq,theta.hat)
    text(x=xmin.Graph,
         y=ymin.Graph+(ymax.Graph-ymin.Graph)*0.90,
         labels=paste("D mean=",D.fixed," cm",sep=""),cex=0.8,pos=4)     
  }
  dev.off()

  ########################################################################################################################
  ########################################################################################################################
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##
  ##                                      Exporting mortality parameters for each species
  ##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ########################################################################################################################
  ########################################################################################################################
  
  th <- th.sign/100 # Threshold for significance
  sppar_mortality <- as.data.frame(matrix(NA,ncol=8,nrow=nspecies+1))
  names(sppar_mortality) <- c("Id","Sp","alpha0m","alpha1m","alpha2m","s0m","s1m","s2m")
  sppar_mortality$Id <- c(0,1:nspecies)
  sppar_mortality$Sp <- c("Spmean",sort(unique(as.character(DataMortality$Sp))))
  #= Loop on species
  for (k in 1:nspecies) {
    sppar_mortality$alpha0m[k+1] <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")])*CF_mortality
    sppar_mortality$alpha1m[k+1] <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])*CF_mortality
    sppar_mortality$alpha2m[k+1] <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])*CF_mortality
    # Significance
    q0 <- quantile(MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_mortality$s0m[k+1] <- ifelse(q0[1]*q0[2]<0,0,1)
    q1 <- quantile(MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_mortality$s1m[k+1] <- ifelse(q1[1]*q1[2]<0,0,1)
    q2 <- quantile(MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_mortality$s2m[k+1] <- ifelse(q2[1]*q2[2]<0,0,1)
  }
  #= Mean species
  sppar_mortality$alpha0m[1] <- mean(MCMC$beta0)*CF_mortality
  sppar_mortality$alpha1m[1] <- mean(MCMC$beta1)*CF_mortality
  sppar_mortality$alpha2m[1] <- mean(MCMC$beta2)*CF_mortality
  # Significance
  q0 <- quantile(MCMC$beta0,c(th/2,1-th/2))
  sppar_mortality$s0m[1] <- ifelse(q0[1]*q0[2]<0,0,1)
  q1 <- quantile(MCMC$beta1,c(th/2,1-th/2))
  sppar_mortality$s1m[1] <- ifelse(q1[1]*q1[2]<0,0,1)
  q2 <- quantile(MCMC$beta2,c(th/2,1-th/2))
  sppar_mortality$s2m[1] <- ifelse(q2[1]*q2[2]<0,0,1)

  write.table(sppar_mortality,file="sppar_mortality.txt",sep="\t",row.names=FALSE)
  
cat("\nFunction entry_mortality has finished, please check files:
1. time_mortality.txt
2. posteriors_mortality.pdf
3. MCMC_mortality.txt
4. theta_dbh.pdf
5. theta_comp.pdf
6. sppar_mortality.txt
in the working directory\n\n") 
flush.console()

}
######################################################################################
