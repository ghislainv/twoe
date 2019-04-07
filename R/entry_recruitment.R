entry_recruitment <- function (data_recruitment,burnin=5000,mcmc=5000,thin=5,th.sign=5) {

  
  # Some necessary objects
  DataRecruitment <- data_recruitment
  DataRecruitment$IdentSpecies <- as.numeric(as.factor(as.character(DataRecruitment$Sp)))
  Levels.Species <- sort(unique(DataRecruitment$IdentSpecies))
  nspecies <- length(Levels.Species)
  Levels.Species.Name <- sort(unique(as.character(DataRecruitment$Sp)))

  # Transformed variables
  BAsp <- DataRecruitment$BAsp_t-0.5
  C <- DataRecruitment$C_t-20

  ###############
  # Gibbs sampler
  
  # Y, Number of recruits on the cell
  # X1, (DataRecruitment$BAsp_t-0.5) # m2.ha-1
  # X2, (DataRecruitment$C_t-20) # m2.ha-1
   
  Start.C.Process <- Sys.time()
  model <-  entry_recruitment_gibbs(fixed=R_tp1~BAsp+C,
                                    random=~BAsp+C,
                                    group="Sp",
                                    interval=DataRecruitment$interval_tp1,
                                    area=DataRecruitment$SCell,
                                    data=DataRecruitment,
                                    burnin=burnin, mcmc=mcmc, thin=thin,
                                    verbose=1, seed=NA,
                                    beta.start=c(0,0,0), sigma2.start=1, Vb.start=diag(1,3), mubeta=c(0,0,0),
                                    Vbeta=diag(1.0E6,3), r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001,
                                    FixOD=1)
  End.C.Process <- Sys.time()
  
  ########################
  # Computation time
  sink(file="time_recruitment.txt")
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
  pdf(file="posteriors_recruitment.pdf")
  plot(mcmc.obj)
  dev.off()

  ########################
  # Reorganizing the MCMCs
  write.table(MCMC,file="MCMC_recruitment.txt",sep="\t",row.names=FALSE)
  # MCMC <- read.table(file="MCMC_recruitment.txt",sep="\t",header=TRUE)
 
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

  # Correction factor
  # See Chave2005, Oecologia and Hadfield's course notes, page 38 and fig. 2.5.
  # http://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf
  # CF <- MCMC$sigma2/2

  # Maximum likelihood approach
  Poisson.Likelihood <- function (lambda,R,Year,SCell) {
    logvrais <- sum(R*log(lambda*Year*SCell)-(lambda*Year*SCell)-log(factorial(R)))
    if (is.na(logvrais)) {logvrais <- -1.0E6}
    else {
      if (logvrais==Inf) {logvrais <- 1.0E6}
      if (logvrais==-Inf) {logvrais <- -1.0E6}
    } 
    return(logvrais)
  }

  # BAsp
  pdf("recruitment_BAsp.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1),las=3)

  for (k in 1:nspecies) {
    
    # C fixed to the mean of the population
    C.fixed <- round(mean(DataRecruitment$C_t[DataRecruitment$IdentSpecies==Levels.Species[k]]))

    # Row data
    BA.sp <- c(DataRecruitment$BAsp_t[DataRecruitment$IdentSpecies==Levels.Species[k]])
    i.sp <- c(DataRecruitment$interval_tp1[DataRecruitment$IdentSpecies==Levels.Species[k]])
    R.sp <- c(DataRecruitment$R_tp1[DataRecruitment$IdentSpecies==Levels.Species[k]])
    S.sp <- c(DataRecruitment$SCell[DataRecruitment$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(BA.sp); xmax.Data <- max(BA.sp)
    
    # Bins for MLL mean number of recruits by Comp class
    bins <- c(0,0.25,0.5,1,5,10,max(20,BA.sp+1))
    midpoints <- c(0.125,0.375,0.75,3,7.5,10+(max(20,BA.sp+1)-10)/2)
    BAsplit <- cut(BA.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    levels.bin <- levels(BAsplit)
    MLL.lambda <- rep(0,length(midpoints))
    for (i in 1:length(midpoints)) {
      R <- R.sp[BAsplit==levels.bin[i]&!is.na(BAsplit)]
      Year <- i.sp[BAsplit==levels.bin[i]&!is.na(BAsplit)]
      SCell <- S.sp[BAsplit==levels.bin[i]&!is.na(BAsplit)]
      MinMLL <- optimize(f=Poisson.Likelihood,lower=1.0E-10,upper=10,R=R,Year=Year,SCell=SCell,maximum=TRUE)
      if (length(Year)==0) {MLL.lambda[i] <- NA}
      if (length(Year)!=0&sum(R)==0) {MLL.lambda[i] <- 0}
      if (length(Year)!=0&sum(R)!=0) {MLL.lambda[i] <- MinMLL$maximum}
    }
    ymin.MLL <- min(MLL.lambda,na.rm=TRUE); ymax.MLL <- max(MLL.lambda,na.rm=TRUE);
    index1 <- which(!is.na(MLL.lambda))[1]
    index2 <- which(!is.na(MLL.lambda))[length(which(!is.na(MLL.lambda)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    BAsp.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)
 
    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    lambda.hat <- vector()
    for (n in 1:100) {
      lambda.hat[n] <- exp(b0sp+b1sp*(BAsp.seq[n]-0.5)+b2sp*(C.fixed-20))
    }
    ymin.Pred <- min(lambda.hat); ymax.Pred <- max(lambda.hat)
 
    # Prediction mean species
    b0.msp <- mean(MCMC$beta0+MCMC$sigma2/2)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    lambda.hat.msp <- vector()
    for (n in 1:100) {
      lambda.hat.msp[n] <- exp(b0.msp+b1.msp*(BAsp.seq[n]-0.5)+b2.msp*(C.fixed-20))
    }
    ymin.Pred.msp <- min(lambda.hat.msp); ymax.Pred.msp <- max(lambda.hat.msp)
    ymin.Graph <- min(ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.MLL,ymax.Pred,ymax.Pred.msp)

    plot(midpoints,MLL.lambda,
         pch=18,
         cex=0.8,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="BA conspecific (m2.ha-1)",
         ylab="Recruits (.yr-1.m-2)")
    lines(BAsp.seq,lambda.hat.msp,col="grey")
    lines(BAsp.seq,lambda.hat)
    text(x=xmin.Graph,
         y=ymin.Graph+(ymax.Graph-ymin.Graph)*0.90,
         labels=paste("C mean=",C.fixed," m2.ha-1",sep=""),cex=0.8,pos=4)
  }
  dev.off()

  
  # C
  pdf("recruitment_comp.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1),las=3)
 
  for (k in 1:nspecies) {
  
    # BAsp fixed to the mean of the population
    BAsp.fixed <- round(mean(DataRecruitment$BAsp_t[DataRecruitment$IdentSpecies==Levels.Species[k]]),2)

    # Row data
    C.sp <- c(DataRecruitment$C_t[DataRecruitment$IdentSpecies==Levels.Species[k]])
    i.sp <- c(DataRecruitment$interval_tp1[DataRecruitment$IdentSpecies==Levels.Species[k]])
    R.sp <- c(DataRecruitment$R_tp1[DataRecruitment$IdentSpecies==Levels.Species[k]])
    S.sp <- c(DataRecruitment$SCell[DataRecruitment$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(C.sp); xmax.Data <- max(C.sp)

    # Bins for MLL mean number of recruits by Comp class
    bins <- c(0,10,20,30,50,70,max(100,C.sp+1))
    midpoints <- c(5,15,25,45,60,70+(max(100,C.sp+1)-70)/2)
    Csplit <-  cut(C.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    levels.bin <- levels(Csplit)
    MLL.lambda <- rep(0,length(midpoints))
    for (i in 1:length(midpoints)) {
      R <- R.sp[Csplit==levels.bin[i]]
      Year <- i.sp[Csplit==levels.bin[i]]
      SCell <- S.sp[Csplit==levels.bin[i]]
      MinMLL <- optimize(f=Poisson.Likelihood,lower=1.0E-10,upper=10,R=R,Year=Year,SCell=SCell,maximum=TRUE)
      if (length(Year)==0) {MLL.lambda[i] <- NA}
      if (length(Year)!=0&sum(R)==0) {MLL.lambda[i] <- 0}
      if (length(Year)!=0&sum(R)!=0) {MLL.lambda[i] <- MinMLL$maximum}
    }
    ymin.MLL <- min(MLL.lambda,na.rm=TRUE); ymax.MLL <- max(MLL.lambda,na.rm=TRUE);
    index1 <- which(!is.na(MLL.lambda))[1]
    index2 <- which(!is.na(MLL.lambda))[length(which(!is.na(MLL.lambda)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    C.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)
 
    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    lambda.hat <- vector()
    for (n in 1:100) {
      lambda.hat[n] <- exp(b0sp+b1sp*(BAsp.fixed-0.5)+b2sp*(C.seq[n]-20))
    }
    ymin.Pred <- min(lambda.hat); ymax.Pred <- max(lambda.hat)
 
    # Prediction mean species
    b0.msp <- mean(MCMC$beta0+MCMC$sigma2/2)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    G.hat.msp <- vector()
    lambda.hat.msp <- vector()
    for (n in 1:100) {
      lambda.hat.msp[n] <- exp(b0.msp+b1.msp*(BAsp.fixed-0.5)+b2.msp*(C.seq[n]-20))
    }
    ymin.Pred.msp <- min(lambda.hat.msp); ymax.Pred.msp <- max(lambda.hat.msp)
    ymin.Graph <- min(ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.MLL,ymax.Pred,ymax.Pred.msp)

    plot(midpoints,MLL.lambda,
         pch=18,
         cex=0.8,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="Competition (m2.ha-1)",
         ylab="Recruits (.yr-1.m-2)")
    lines(C.seq,lambda.hat.msp,col="grey")
    lines(C.seq,lambda.hat)
    text(x=xmin.Graph,
         y=ymin.Graph+(ymax.Graph-ymin.Graph)*0.90,
         labels=paste("BAsp mean=",BAsp.fixed," m2.ha-1",sep=""),cex=0.8,pos=4)
  }
  dev.off()

  ########################################################################################################################
  ########################################################################################################################
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##
  ##                                      Exporting recruitment parameters for each species
  ##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ########################################################################################################################
  ########################################################################################################################
  
  th <- th.sign/100 # Threshold for significance
  sppar_recruitment <- as.data.frame(matrix(NA,ncol=8,nrow=nspecies+1))
  names(sppar_recruitment) <- c("Id","Sp","alpha0r","alpha1r","alpha2r","s0r","s1r","s2r")
  sppar_recruitment$Id <- c(0,1:nspecies)
  sppar_recruitment$Sp <- c("Spmean",sort(unique(as.character(DataRecruitment$Sp))))
  #= Loop on species
  for (k in 1:nspecies) {
    sppar_recruitment$alpha0r[k+1] <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    sppar_recruitment$alpha1r[k+1] <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    sppar_recruitment$alpha2r[k+1] <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    # Significance
    q0 <- quantile(MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_recruitment$s0r[k+1] <- ifelse(q0[1]*q0[2]<0,0,1)
    q1 <- quantile(MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_recruitment$s1r[k+1] <- ifelse(q1[1]*q1[2]<0,0,1)
    q2 <- quantile(MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_recruitment$s2r[k+1] <- ifelse(q2[1]*q2[2]<0,0,1)
  }
  sppar_recruitment$alpha0r[1] <- mean(MCMC$beta0+MCMC$sigma2/2)
  sppar_recruitment$alpha1r[1] <- mean(MCMC$beta1)
  sppar_recruitment$alpha2r[1] <- mean(MCMC$beta2)
  # Significance
  q0 <- quantile(MCMC$beta0,c(th/2,1-th/2))
  sppar_recruitment$s0r[1] <- ifelse(q0[1]*q0[2]<0,0,1)
  q1 <- quantile(MCMC$beta1,c(th/2,1-th/2))
  sppar_recruitment$s1r[1] <- ifelse(q1[1]*q1[2]<0,0,1)
  q2 <- quantile(MCMC$beta2,c(th/2,1-th/2))
  sppar_recruitment$s2r[1] <- ifelse(q2[1]*q2[2]<0,0,1)
  
  write.table(sppar_recruitment,file="sppar_recruitment.txt",sep="\t",row.names=FALSE)
  
cat("\nFunction entry_recruitment has finished, please check files:
1. time_recruitment.txt
2. posteriors_recruitment.pdf
3. MCMC_recruitment.txt
4. recruitment_dbh.pdf
5. recruitment_comp.pdf
6. sppar_recruitment.txt
in the working directory\n\n") 
flush.console()
  
}
######################################################################################
