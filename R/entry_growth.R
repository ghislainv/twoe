entry_growth <- function (data_growth,burnin=1000,mcmc=1000,thin=1,th.sign=5) {

  # Some necessary objects
  DataGrowth <- data_growth[data_growth$G_tp1>(-2)&data_growth$G_tp1<100,]
  DataGrowth$IdentSpecies <- as.numeric(as.factor(as.character(DataGrowth$Sp)))
  Levels.Species <- sort(unique(DataGrowth$IdentSpecies))
  nspecies <- length(Levels.Species)
  Levels.Species.Name <- sort(unique(as.character(DataGrowth$Sp)))

  # Transformed variables
  lG <- log(DataGrowth$G_tp1+2)
  lD <- log(DataGrowth$D_t)
  lC <- log(DataGrowth$C_t+1)

  ###############
  # Gibbs sampler
  
  # Y, Diameter growth in mm.yr-1
  # X1, log(DataGrowth$D_t) # cm
  # X2, log(DataGrowth$C_t +1) # m2.ha-1
  
  Start.C.Process <- Sys.time()
  model <- entry_growth_gibbs(fixed=lG~lD+lC,
                              random=~lD+lC,
                              group="Sp",
                              data=DataGrowth,
                              burnin=burnin, mcmc=mcmc, thin=thin,
                              verbose=1, seed=NA,
                              beta.start=c(0,0,0), sigma2.start=1, Vb.start=diag(1,3), mubeta=c(0,0,0),
                              Vbeta=diag(1.0E6,3), r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001)
  End.C.Process <- Sys.time()
 
  ########################
  # Computation time
  sink(file="time_growth.txt")
  Time.C.Process <- difftime(End.C.Process,Start.C.Process,units="mins")
  cat("N# iterations: ",burnin+mcmc,"\n","Duration of the C++ program (in mins): ",Time.C.Process,sep="")
  sink()
 
  ##################################################################################################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #                                     Analysing MCMC ouputs                                      #
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##################################################################################################

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
  pdf(file="posteriors_growth.pdf")
  plot(mcmc.obj)
  dev.off()
  
  ########################
  # Export MCMCs
  write.table(MCMC,file="MCMC_growth.txt",sep="\t",row.names=FALSE)
  # MCMC <- read.table(file="MCMC_growth.txt",sep="\t",header=TRUE)
 
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
  
  # D
  pdf("growth_dbh.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))
  
  for (k in 1:nspecies) {
    
    # C fixed to the mean of the population
    C.fixed <- round(mean(DataGrowth$C_t[DataGrowth$IdentSpecies==Levels.Species[k]]))

    # Row data
    D.sp <- c(DataGrowth$D_t[DataGrowth$IdentSpecies==Levels.Species[k]])
    G.sp <- c(DataGrowth$G_tp1[DataGrowth$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(D.sp); xmax.Data <- max(D.sp)
    ymin.Data <- min(G.sp); ymax.Data <- max(G.sp)

    # Bins for arithmetic mean
    bins <- c(min(10,D.sp-1),20,40,60,100,max(200,D.sp+1))
    midpoints <- c(min(10,D.sp-1)+(20-min(10,D.sp-1))/2,30,50,80,100+(max(200,D.sp+1)-100)/2)
    dbhsplit <- cut(D.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    dincmeans <- as.numeric(tapply(G.sp,dbhsplit,mean,na.rm=TRUE))
    ymin.MLL <- min(dincmeans,na.rm=TRUE); ymax.MLL <- max(dincmeans,na.rm=TRUE);
    index1 <- which(!is.na(dincmeans))[1]
    index2 <- which(!is.na(dincmeans))[length(which(!is.na(dincmeans)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    D.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)
 
    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    G.hat <- vector()
    for (n in 1:100) {
      G.hat[n] <- exp(b0sp+b1sp*log(D.seq[n])+b2sp*log(C.fixed+1))-2
    }
    ymin.Pred <- min(G.hat); ymax.Pred <- max(G.hat)

    # Prediction mean species
    b0.msp <- mean(MCMC$beta0+MCMC$sigma2/2)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    G.hat.msp <- vector()
    for (n in 1:100) {
      G.hat.msp[n] <- exp(b0.msp+b1.msp*log(D.seq[n])+b2.msp*log(C.fixed+1))-2
    }
    ymin.Pred.msp <- min(G.hat.msp); ymax.Pred.msp <- max(G.hat.msp)
    ymin.Graph <- min(ymin.Data,ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.Data,ymax.MLL,ymax.Pred,ymax.Pred.msp)

    plot(D.sp,G.sp,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="DBH (cm)",
         ylab="DBH growth (mm.yr-1)",
         col="grey",
         cex=0.6)
    points(midpoints,dincmeans,col="black", pch=18, cex=0.8)
    lines(D.seq,G.hat.msp,col="grey")
    lines(D.seq,G.hat)
    text(x=xmin.Graph,
         y=ymin.Graph+(ymax.Graph-ymin.Graph)*0.90,
         labels=paste("C mean=",C.fixed," m2.ha-1",sep=""),cex=0.8,pos=4)
  }
  dev.off()
 
  # C
  pdf("growth_comp.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(4,3),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))
 
  for (k in 1:nspecies) {

    # D fixed to the mean of the population
    D.fixed <- round(mean(DataGrowth$D_t[DataGrowth$IdentSpecies==Levels.Species[k]]))
    
    # Row data
    C.sp <- c(DataGrowth$C_t[DataGrowth$IdentSpecies==Levels.Species[k]])
    G.sp <- c(DataGrowth$G_tp1[DataGrowth$IdentSpecies==Levels.Species[k]])
    xmin.Data <- min(C.sp); xmax.Data <- max(C.sp)
    ymin.Data <- min(G.sp); ymax.Data <- max(G.sp)

    # Bins for arithmetic mean
    bins <- c(0,10,20,30,50,70,max(100,C.sp+1))
    midpoints <- c(5,15,25,45,60,70+(max(100,C.sp+1)-70)/2)
    Csplit <-  cut(C.sp,breaks=bins,right=FALSE,include.lowest=TRUE)
    dincmeans <- as.numeric(tapply(G.sp,Csplit,mean,na.rm=TRUE))
    index1 <- which(!is.na(dincmeans))[1]
    index2 <- which(!is.na(dincmeans))[length(which(!is.na(dincmeans)))]
    xmin.MLL <- min(midpoints[index1]); xmax.MLL <- max(midpoints[index2])
    xmin.Graph <- min(xmin.Data,xmin.MLL); xmax.Graph <- max(xmax.Data,xmax.MLL)
    C.seq <- seq(from=xmin.Graph,to=xmax.Graph,length.out=100)
 
    # Prediction
    b0sp <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    b1sp <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    b2sp <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    G.hat <- vector()
    for (n in 1:100) {
      G.hat[n] <- exp(b0sp+b1sp*log(D.fixed)+b2sp*log(C.seq[n]+1))-2
    }
    ymin.Pred <- min(G.hat); ymax.Pred <- max(G.hat)

    # Prediction mean species
    b0.msp <- mean(MCMC$beta0+MCMC$sigma2/2)
    b1.msp <- mean(MCMC$beta1)
    b2.msp <- mean(MCMC$beta2)
    G.hat.msp <- vector()
    for (n in 1:100) {
      G.hat.msp[n] <- exp(b0.msp+b1.msp*log(D.fixed)+b2.msp*log(C.seq[n]+1))-2
    }
    ymin.Pred.msp <- min(G.hat.msp); ymax.Pred.msp <- max(G.hat.msp)
    ymin.Graph <- min(ymin.Data,ymin.MLL,ymin.Pred,ymin.Pred.msp)
    ymax.Graph <- max(ymax.Data,ymax.MLL,ymax.Pred,ymax.Pred.msp)
    
    plot(C.sp,G.sp,
         xlim=c(xmin.Graph,xmax.Graph),ylim=c(ymin.Graph,ymax.Graph),
         main=Levels.Species.Name[k],
         xlab="Competition (m2.ha-1)",
         ylab="DBH growth (mm.yr-1)",
         col="grey",
         cex=0.6)
    points(midpoints,dincmeans,col="black", pch=18, cex=0.8)
    lines(C.seq,G.hat.msp,col="grey")
    lines(C.seq,G.hat)
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
  ##                                      Exporting growth parameters for each species
  ##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
  ########################################################################################################################
  ########################################################################################################################

  th <- th.sign/100 # Threshold for significance
  sppar_growth <- as.data.frame(matrix(NA,ncol=8,nrow=nspecies+1))
  names(sppar_growth) <- c("Id","Sp","alpha0g","alpha1g","alpha2g","s0g","s1g","s2g")
  sppar_growth$Id <- c(0,1:nspecies)
  sppar_growth$Sp <- c("Spmean",sort(unique(as.character(DataGrowth$Sp))))
  #= Loop on species
  for (k in 1:nspecies) {
    sppar_growth$alpha0g[k+1] <- mean(MCMC$beta0+MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")]+MCMC$sigma2/2)
    sppar_growth$alpha1g[k+1] <- mean(MCMC$beta1+MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")])
    sppar_growth$alpha2g[k+1] <- mean(MCMC$beta2+MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")])
    # Significance
    q0 <- quantile(MCMC[,names(MCMC)==paste("b0",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_growth$s0g[k+1] <- ifelse(q0[1]*q0[2]<0,0,1)
    q1 <- quantile(MCMC[,names(MCMC)==paste("b1",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_growth$s1g[k+1] <- ifelse(q1[1]*q1[2]<0,0,1)
    q2 <- quantile(MCMC[,names(MCMC)==paste("b2",Levels.Species[k],sep=".")],c(th/2,1-th/2))
    sppar_growth$s2g[k+1] <- ifelse(q2[1]*q2[2]<0,0,1)
  }
  #= Mean species
  sppar_growth$alpha0g[1] <- mean(MCMC$beta0+MCMC$sigma2/2)
  sppar_growth$alpha1g[1] <- mean(MCMC$beta1)
  sppar_growth$alpha2g[1] <- mean(MCMC$beta2)
  # Significance
  q0 <- quantile(MCMC$beta0,c(th/2,1-th/2))
  sppar_growth$s0g[1] <- ifelse(q0[1]*q0[2]<0,0,1)
  q1 <- quantile(MCMC$beta1,c(th/2,1-th/2))
  sppar_growth$s1g[1] <- ifelse(q1[1]*q1[2]<0,0,1)
  q2 <- quantile(MCMC$beta2,c(th/2,1-th/2))
  sppar_growth$s2g[1] <- ifelse(q2[1]*q2[2]<0,0,1)
  
  write.table(sppar_growth,file="sppar_growth.txt",sep="\t",row.names=FALSE)
  
cat("\nFunction entry_growth has finished, please check files:
1. time_growth.txt
2. posteriors_growth.pdf
3. MCMC_growth.txt
4. growth_dbh.pdf
5. growth_comp.pdf
6. sppar_growth.txt
in the working directory\n\n") 
flush.console()
  
}
######################################################################################
