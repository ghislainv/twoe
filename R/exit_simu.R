exit_simu <- function(Data,XY.Plot,D.Recruitment,R.Comp,L.Cell,Plot.Sim,Year.Sim,Step.Sim) {

  #============================================================#
  # Time for each census (including 0)
  #============================================================#
  nCens <- ncol(Data)-5 # Number of censuses 
  Date.dash <- rep(NA,nCens)
  Interval <- rep(NA,nCens-1)
  for (d in 1:nCens) {
    date.point <- strsplit(x=names(Data)[5+d],split="D")[[1]][2]
    Date.dash[d] <- gsub("\\.","-",date.point)
  }
  for (d in 1:(nCens-1)) {
    Interval[d] <- difftime(Date.dash[d+1],Date.dash[d],units="days")/365.25
  }
  Time <- c(0,Interval)
  for (d in 2:(nCens)) {
     Time[d] <- Time[d]+Time[d-1]
  }
  
  #============================================================#
  # Merging data-sets of species parameters
  #============================================================#

  sppar_growth <- read.table(file="sppar_growth.txt",header=TRUE,sep="\t")
  sppar_mortality <- read.table(file="sppar_mortality.txt",header=TRUE,sep="\t")
  sppar_recruitment <- read.table(file="sppar_recruitment.txt",header=TRUE,sep="\t")
  sppar_gm <- merge(sppar_growth,sppar_mortality,by="Sp",all=TRUE,sort=TRUE)
  sppar <- merge(sppar_gm,sppar_recruitment,by="Sp",all=TRUE,sort=TRUE)

  #== Some species may not have parameter depending on the process
  #== We affect them the parameters of the mean species
  for (i in 1:ncol(sppar)) {
    sppar[is.na(sppar[,i]),i] <- sppar[sppar$Sp=="Spmean",i]
  }

  #== Dealing with non-significant parameters
  for (i in 1:nrow(sppar)) {
    sppar[i,c(3:5)] <- sppar[i,c(3:5)]*sppar[i,c(6:8)]
    sppar[i,c(10:12)] <- sppar[i,c(10:12)]*sppar[i,c(13:15)]
    sppar[i,c(17:19)] <- sppar[i,c(17:19)]*sppar[i,c(20:22)]
  }
  for (i in c(3:5,10:12,17:19)) {
    sppar[sppar[,i]==0,i] <- sppar[sppar$Sp=="Spmean",i]
  }
 
  #============================================================#
  # Selecting the plots for simulations
  #============================================================#

  #===== Number of plots
  nPlot.Sim <- length(Plot.Sim)
  
  #===== Loop on plots
  for (j in 1:nPlot.Sim) {

    #===== Dimensions of the plot
    X.Plot <- XY.Plot$X.Plot[XY.Plot$Plot==Plot.Sim[j]]
    Y.Plot <- XY.Plot$Y.Plot[XY.Plot$Plot==Plot.Sim[j]]
    
    #===== Message
    if (j==1) {
      cat(paste("\nProcessing plot ",Plot.Sim[j],"...\n",sep=""))
      flush.console()
    }
    if (j!=1) {
      cat(paste("Processing plot ",Plot.Sim[j],"...\n",sep=""))
      flush.console()	
    }

    #===== Selecting initial data for the plot
    data_sim <- Data[Data$Plot==Plot.Sim[j],c(1:6)]
    data_sim[data_sim==-9999] <- 0
    names(data_sim)[6] <- "D0"
    
    #===== Number of cells on each plot depending on plot size (X.Plot,Y.Plot)
    nXCell <- floor(X.Plot/L.Cell)
    nYCell <- floor(Y.Plot/L.Cell)
    nCell <- nXCell*nYCell
    Lx.LastCell <- L.Cell+X.Plot%%L.Cell
    Ly.LastCell <- L.Cell+Y.Plot%%L.Cell

    #===== Data.Cell for recruitment process
    Data.Cell <- data.frame(matrix(NA,ncol=7,nrow=nCell))
    colnames(Data.Cell) <- c("Plot","Cell","X.Cell.min","X.Cell.max","Y.Cell.min","Y.Cell.max","S.Cell")
    Data.Cell$Plot <- Plot.Sim[j]
    Data.Cell$Cell <- c(1:nCell)
    X.Cell.min <- c(seq(from=0,to=(nXCell-1)*L.Cell,by=L.Cell))
    X.Cell.max <- c(seq(from=L.Cell,to=(nXCell-2)*L.Cell+L.Cell,by=L.Cell),X.Plot)
    Y.Cell.min <- c(seq(from=0,to=(nYCell-1)*L.Cell,by=L.Cell))
    Y.Cell.max <- c(seq(from=L.Cell,to=(nYCell-2)*L.Cell+L.Cell,by=L.Cell),Y.Plot)
    Data.Cell$X.Cell.min <- rep(X.Cell.min,nYCell)
    Data.Cell$X.Cell.max <- rep(X.Cell.max,nYCell)
    Data.Cell$Y.Cell.min <- rep(Y.Cell.min,each=nXCell)
    Data.Cell$Y.Cell.max <- rep(Y.Cell.max,each=nXCell)
    Data.Cell$S.Cell <- (Data.Cell$X.Cell.max-Data.Cell$X.Cell.min)*(Data.Cell$Y.Cell.max-Data.Cell$Y.Cell.min) # S in m2

    #===== Initial number for new tree (recruitment process)
    twoetree <- 0
    
    #===== Loop on years of simulations
    for (d in 1:Year.Sim) {

      #========================================================#
      # Matrix of parameters
      #========================================================#
      # Be careful with the use of merge and the ordering !!
      data_sim$order <- seq(len=nrow(data_sim))
      m <- merge(sppar,data_sim,by.x="Sp",by.y="Species",all.y=TRUE,sorted=FALSE)
      MofPar <- m[sort.list(m$order),c(3:5,10:12,17:19)]
      data_sim <- data_sim[,-c(ncol(data_sim))]
      
      #========================================================#
      # Computing the competition index for each tree
      #========================================================#

      #===== Number of trees for plot
      Levels.IdentTree <- data_sim$Tree
      nTree.Plot <- length(data_sim$Tree)

      #===== Computing the distance between trees with dist function
      Dist.Matrix <- as.matrix(dist(cbind(data_sim$X,data_sim$Y),diag=TRUE,upper=TRUE))

      #===== Trees with a distance <= to R.Comp
      Include.Mat <- Dist.Matrix
      Include.Mat[Include.Mat<=R.Comp] <- 1
      Include.Mat[Include.Mat>R.Comp] <- 0

      #===== We do not include the target tree: 0 on diagonal
      for (i in 1:nTree.Plot) {
        Include.Mat[i,i] <- 0
      }
  
      #===== Basal Area vectors
      D.vect <- data_sim[,ncol(data_sim)]
      BA.vect <- pi*(D.vect/200)^2

      #===== Surface of the neighborhood for each tree
      Sm2 <- surface_comp_index(X.trees=data_sim$X,
                                Y.trees=data_sim$Y,
                                X.Plot=X.Plot,
                                Y.Plot=Y.Plot,
                                R.Comp=R.Comp)
      Sha <- Sm2/10000 # in ha

      #===== Competition index for each tree
      CI.vect <- Include.Mat%*%BA.vect/Sha

      #========================================================#
      # Growth process
      #========================================================#     
      
      G.vect <- exp(MofPar$alpha0g+MofPar$alpha1g*log(D.vect)+MofPar$alpha2g*log(CI.vect+1))-2
      eval(parse(text=paste("data_sim$D",d," <- 0",sep="")))
      eval(parse(text=paste("data_sim$D",d," <- data_sim$D",d-1,"+as.numeric(G.vect)/10",sep="")))
      data_sim[data_sim[ncol(data_sim)-1]==0,ncol(data_sim)] <- 0 
      
      #========================================================#
      # Mortality process
      #========================================================#

      #===== Probability of dying
      theta.vect <-  inv_logit(MofPar$alpha0m+MofPar$alpha1m*(D.vect-20)+MofPar$alpha2m*(CI.vect-20))
      z <- runif(nTree.Plot,0,1)
      # z[D.vect>250]==0 # We assume that trees with D>250cm must die

      #===== Set D=0 to dead trees
      data_sim[z<=theta.vect,ncol(data_sim)] <- 0

      #========================================================#
      # Recruitment process
      #========================================================#

      #===== Number of species
      Levels.Species <- sort(unique(as.factor(as.character(data_sim$Sp))))
      nSpecies.Plot <- length(Levels.Species)
      IdentSpecies <- as.numeric(as.factor(as.character(data_sim$Sp)))

      #===== Trees on each cell
      TreesOnCell <- vector()
      for (i in 1:nTree.Plot) {
        TreesOnCell[i] <- min(which(data_sim$X[i]>=Data.Cell$X.Cell.min&
                                    data_sim$X[i]<=Data.Cell$X.Cell.max&
                                    data_sim$Y[i]>=Data.Cell$Y.Cell.min&
                                    data_sim$Y[i]<=Data.Cell$Y.Cell.max)) # min allocates the tree in a frontier to the
                                                                           # lower cell number
      }

      #===== Competition index on each cell
      BAr.vect <- pi*(D.vect/200)^2
      for (c in 1:nCell) {
         Data.Cell$C[c] <- c(as.integer(TreesOnCell==c)%*%BAr.vect)/(Data.Cell$S.Cell[c]/10000) # Surface in ha
      }

      #===== BAsp and lambda on each cell
      BAsp <- lambda <- matrix(NA,nrow=nCell,ncol=nSpecies.Plot)
      colnames(BAsp) <- colnames(lambda) <- Levels.Species
      for (c in 1:nCell) {
        for (k in 1:nSpecies.Plot) {
          #== BAsp
          BAsp[c,colnames(BAsp)==Levels.Species[k]] <-
            c(as.integer(TreesOnCell==c&data_sim$Sp==as.character(Levels.Species[k]))%*%BAr.vect)/(Data.Cell$S.Cell[c]/10000)
          #== lambda
          alpha0r <- sppar_recruitment$alpha0r[sppar_recruitment$Sp==as.character(Levels.Species[k])]
          alpha1r <- sppar_recruitment$alpha1r[sppar_recruitment$Sp==as.character(Levels.Species[k])]
          alpha2r <- sppar_recruitment$alpha2r[sppar_recruitment$Sp==as.character(Levels.Species[k])]
          lambda[c,colnames(lambda)==Levels.Species[k]] <- exp(alpha0r+
                             alpha1r*(BAsp[c,colnames(BAsp)==Levels.Species[k]]-0.5)+
                             alpha2r*(Data.Cell$C[c]-20))*Data.Cell$S.Cell[c] # must be multiplied by the surface of the cell in m2
        }
      }
      
      #===== Rsp
      Rsp <- matrix(rpois(nCell*nSpecies.Plot,c(lambda)),nrow=nCell,ncol=nSpecies.Plot)
      colnames(Rsp) <- Levels.Species

      #===== Constructing a new data for recruits
      data_newtree <-  data_sim[0,]
      for (c in 1:nCell) {
        for (k in 1:nSpecies.Plot) {
          nRecruit <- as.numeric(Rsp[c,k])
          #== Only if there are recruits on cell c for species k
          if (nRecruit>0) {
            #== newtree matrix
            newtree <- data.frame(matrix(NA,ncol=ncol(data_sim),nrow=nRecruit))
            names(newtree) <- names(data_sim)
            #== We fill in newtree matrix
            newtree$Species <- Levels.Species[k]
            newtree$Plot <- rep(Plot.Sim[j],nRecruit)
            #== New tree id from twoetree index
            seqnewtree <- c((twoetree+1):(twoetree+nRecruit))
            twoetree <- twoetree+nRecruit
            newtree$Tree <- paste("twoe",seqnewtree,sep="")
            #== newtree coordinates
            newtree$X <- runif(nRecruit,Data.Cell$X.Cell.min[c],Data.Cell$X.Cell.max[c])
            newtree$Y <- runif(nRecruit,Data.Cell$Y.Cell.min[c],Data.Cell$Y.Cell.max[c])
            #== diameters values
            newtree[,6:(ncol(data_sim)-1)] <- 0
            newtree[,ncol(data_sim)] <- D.Recruitment
            #== adding these recruits to data_newtree
            data_newtree <- rbind(data_newtree,newtree)
          }
        }
      }

      #===== Adding new recruits at each iterative year to data_simu
      data_sim <- rbind(data_sim,data_newtree)

      #===== We keep only the first, the two last columns and years each Step.Sim
      lds <- ncol(data_sim)
      if (d>=3 & (d-2)%%Step.Sim!=0) {
        data_sim <- data_sim[,-c(lds-2)]
      }
      
      #===== Message for year iteration
      if (d==1) {
        cat(paste("Year ",d,",",sep=""))
        flush.console()
      }
      if (d!=1 & d!=Year.Sim) {
        if (d%%20!=0) {
          cat(paste(d,",",sep=""))
          flush.console()	
        }
        if (d%%20==0) {
          cat(paste(d,"\nYear ",sep=""))
          flush.console()	
        }
      }
      if (d==Year.Sim) {
        cat(paste(d,"\n\n",sep=""))
        flush.console()
      }  
    }

    #========================================================#
    # Exporting the data-set including simulations
    #========================================================#

    #===== We do not keep the next to last simulation
    data_sim <- data_sim[,-c(ncol(data_sim)-1)]

    #===== Saving the simulations for each plot
    eval(parse(text=paste("write.table(data_sim,file=\"data_sim_plot_",Plot.Sim[j],".txt\",sep=\"\t\",row.names=FALSE)",sep="")))

  }
 
  #========================================================#
  # Graphics of the evolution of the Basal Area
  #========================================================#
  
  pdf("exit_BA.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(2,2),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))

  #===== Loop on plots
  for (j in 1:nPlot.Sim) {

    #===== Dimensions of the plot
    X.Plot <- XY.Plot$X.Plot[XY.Plot$Plot==Plot.Sim[j]]
    Y.Plot <- XY.Plot$Y.Plot[XY.Plot$Plot==Plot.Sim[j]]
    
    #== Original data
    Data.orig <- Data[Data$Plot==Plot.Sim[j],]
    Data.orig[Data.orig==-9999] <- 0
    BA.orig <- apply(pi*(Data.orig[,6:ncol(Data.orig)]/200)^2,2,sum)/((X.Plot*Y.Plot)/10000)
    y.min.orig <- min(BA.orig); y.max.orig <- max(BA.orig)
    
    #== Importing data_sim
    eval(parse(text=paste("data_sim <- read.table(file=\"data_sim_plot_",Plot.Sim[j],".txt\",header=TRUE,sep=\"\t\")",sep="")))
    #== Computation of the basal area
    BA.plot <- apply(pi*(data_sim[,6:ncol(data_sim)]/200)^2,2,sum)/((X.Plot*Y.Plot)/10000)
    y.min.sim <- min(BA.plot); y.max.sim <- max(BA.plot)
    y.min.Graph <- floor(min(y.min.orig,y.min.sim,10))
    y.max.Graph <- ceiling(max(y.max.orig,y.max.sim,40))
    Time.Sim <- as.numeric(sub("D","",names(BA.plot)))
    
    #== Plot
    plot(x=Time.Sim,y=BA.plot,
         ylim=c(y.min.Graph,y.max.Graph),
         type="l",
         main=paste("Plot ",Plot.Sim[j],sep=""),
         xlab="Time (years)",
         ylab="BA (m2.ha-1)",
         axes=FALSE)
    axis(1,at=Time.Sim,labels=Time.Sim)
    axis(2,at=seq(from=y.min.Graph,to=y.max.Graph,by=5),labels=seq(from=y.min.Graph,to=y.max.Graph,by=5))
    #== Original data
    points(Time,BA.orig,col="black")
    mtext(Date.dash[1], side=1, line=+2, at=Time[1],cex=0.6)
    
  }
  dev.off()

  #============================================================#
  # Graphics of the evolution of the BA by species
  #============================================================#
  
  pdf("exit_biodiversity.pdf", pointsize=8, onefile=TRUE, family="Helvetica", paper="letter")
  par(mfrow=c(2,2),cex=1, las=1, lwd=1, pch=18, cex.lab=1, mar=c(4,4,2,1))

  #===== Loop on plots
  for (j in 1:nPlot.Sim) {

    #===== Dimensions of the plot
    X.Plot <- XY.Plot$X.Plot[XY.Plot$Plot==Plot.Sim[j]]
    Y.Plot <- XY.Plot$Y.Plot[XY.Plot$Plot==Plot.Sim[j]]
    
    #== Original data
    Data.orig <- Data[Data$Plot==Plot.Sim[j],]
    Data.orig[Data.orig==-9999] <- 0

    #== Number of species at date t=0
    Levels.Species <- unique(as.character(Data.orig$Species))
    nSpecies <- length(Levels.Species)
    BA.Sp.orig <- matrix(0,ncol=nCens,nrow=nSpecies)
    for (k in 1:nSpecies) {
      Data.Sp <- Data.orig[Data.orig$Species==Levels.Species[k],6:ncol(Data.orig)]
      BA.Sp.orig[k,] <- apply(pi*(Data.Sp/200)^2,2,sum)/((X.Plot*Y.Plot)/10000)
    }
    y.min.orig <- min(BA.Sp.orig); y.max.orig <- max(BA.Sp.orig)
    
    #== Importing data_sim
    eval(parse(text=paste("data_sim <- read.table(file=\"data_sim_plot_",Plot.Sim[j],".txt\",header=TRUE,sep=\"\t\")",sep="")))
    #== Computation of the basal area
    BA.Sp.sim <- matrix(0,ncol=length(6:ncol(data_sim)),nrow=nSpecies)
    for (k in 1:nSpecies) {
      Data.Sp.sim <- data_sim[data_sim$Species==Levels.Species[k],6:ncol(data_sim)]
      BA.Sp.sim[k,] <- apply(pi*(Data.Sp.sim/200)^2,2,sum)/((X.Plot*Y.Plot)/10000)
    }
    y.min.sim <- min(BA.Sp.sim); y.max.sim <- max(BA.Sp.sim)
    y.min.Graph <- 0
    y.max.Graph <- ceiling(max(y.max.orig,y.max.sim,4))
    Time.Sim <- as.numeric(sub("D","",names(BA.plot)))
    
    #== Plot
    #== Colors
    Colors <- rep("black",nSpecies)
    FastSpecies <- order(BA.Sp.orig[,1]-BA.Sp.orig[,nCens])
    Colors[FastSpecies[1:5]] <- "red"
    #== Names fastest species
    Names.Fast <- Levels.Species[FastSpecies[1:5]]
    #== First species
    plot(x=Time.Sim,y=BA.Sp.sim[1,],
         ylim=c(y.min.Graph,y.max.Graph),
         type="l",
         main=paste("Plot ",Plot.Sim[j],sep=""),
         xlab="Time (years)",
         ylab="BA (m2.ha-1)",
         col=Colors[1],
         axes=FALSE)
    axis(1,at=Time.Sim,labels=Time.Sim)
    axis(2,at=seq(from=y.min.Graph,to=y.max.Graph,by=1),labels=seq(from=y.min.Graph,to=y.max.Graph,by=1))
    #== Original data
    points(Time,BA.Sp.orig[1,],col=Colors[1])
    mtext(Date.dash[1], side=1, line=+2, at=Time[1],cex=0.6)
    #== Other species
    for (k in 1:nSpecies) {
      lines(x=Time.Sim,y=BA.Sp.sim[k,],col=Colors[k])
      points(Time,BA.Sp.orig[k,],col=Colors[k])
    }
    #== Text fast species
    text(x=Time.Sim[1],y=0.875*y.max.Graph,
         labels=paste("Fast growing species:\n\n",Names.Fast[1],", ",Names.Fast[2],"\n",Names.Fast[3],", ",Names.Fast[4],"\n",Names.Fast[5],sep=""),
         adj=0,pos=4,cex=0.8)
  }
  dev.off()

  #========================================================#
  # Message
  #========================================================#
  
  cat("\nFunction exit_simu has finished, please check files:
  1. data_simu_plot_\"plotname\".txt
  2. exit_BA.pdf
  3. exit_biodiversity.pdf
  in the working directory\n\n")
  flush.console()
    
}

