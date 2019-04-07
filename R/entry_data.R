## Notes:
## diameter for dead trees must be set to 0
## diameter for trees which have not recruited yet must be set to 0
## not available diameter must be set to -9999

entry_data <-
function(Data,XY.Plot,R.Comp,L.Cell) {

cat("\nFunction entry_data is running.\n")
flush.console()

# Functions
Last.False.True.f <- function(C.Vector) { # Function which return the index of TRUE for the last FALSE-TRUE pair
  null.Vector <- c(C.Vector==0)
  null.FT.Vector <- vector()
  for (i in 1:(length(null.Vector)-1)){
    null.FT.Vector[i] <- (null.Vector[i]==FALSE & null.Vector[i+1]==TRUE)
  }
  Position <- which(null.FT.Vector==TRUE)
  if (sum(null.FT.Vector)!=0) return(Position[length(Position)]+1) else return(0)
}

First.True.f <- function(C.Vector) { # Function which return the index of the first TRUE 
  null.Vector <- c(C.Vector!=0)
  Position <- which(null.Vector==TRUE)
  return(Position[1])
}

# Some necessary variables
nCens <- dim(Data)[2]-5 # Number of censuses
nTree <- length(Data$Tree) # Number of trees
List.Plot <- sort(unique(Data$Plot)) # List of plots
nPlot <- length(List.Plot) # Number of plots

# Time interval
Date.dash <- rep(NA,nCens)
Interval <- rep(NA,nCens-1)
for (d in 1:nCens) {
  date.point <- strsplit(x=names(Data)[5+d],split="D")[[1]][2]
  Date.dash[d] <- gsub("\\.","-",date.point)
}
for (d in 1:(nCens-1)) {
  Interval[d] <- difftime(Date.dash[d+1],Date.dash[d],units="days")/365.25
}

# "Apply" to obtain the date of death and date of recruitment for each tree
ListD <- Data[,6:(6+nCens-1)]
Data$CensusDeath <- apply(ListD,1,Last.False.True.f)
Data$CensusRecruit <- apply(ListD,1,First.True.f)

# Creation of the object Data.pretreat that will be filled in with all plot data
Data.pretreat.g.m <- data.frame() # For growth and mortality (object=the tree)
Data.pretreat.r <- data.frame() # For recruitment (object=the ground cell)

#====================================#
# Simpler recruitment data
#====================================#

# If the time-interval between two censuses is < 5 yr,
# we merge the censuses for the recruitment observation number

# Which censuses do we keep ? => DateR
IntervalR <- Interval
nCensR <- length(IntervalR)
DateR <- c(2:(nCensR+1))
while (sum(IntervalR<5)>1) {
  for (i in 1:(nCensR-1)) {
    if (IntervalR[i]<5) {
      IntervalR[i+1] <- IntervalR[i]+IntervalR[i+1]
      IntervalR <- IntervalR[-c(i)]
      DateR <-  DateR[-c(i)]
      nCensR <- nCensR-1
      break
    }
  }
}

# List of years to merge
yr.List <- list()
yr.List[[1]] <- c(2:DateR[1])
if (length(DateR)>=2) {
  for (i in 2:length(DateR)) {
    yr.List[[i]] <- c((DateR[i-1]+1):DateR[i])
  }
}


# Loop on the plots
for (j in 1:nPlot) {

  #====================================#
  # Message
  #====================================#
  
  if (j==1) {
    cat(paste("\nProcessing plot ",List.Plot[j],"...\n",sep=""))
    flush.console()
  }
  if (j!=1) {
    cat(paste("Processing plot ",List.Plot[j],"...\n",sep=""))
    flush.console()	
  }

  #====================================#
  # Dimensions of the plot
  #====================================#

  X.Plot <- XY.Plot$X.Plot[XY.Plot$Plot==as.character(List.Plot[j])]
  Y.Plot <- XY.Plot$Y.Plot[XY.Plot$Plot==as.character(List.Plot[j])]

  #====================================#
  # Preparing Growth and Mortality data
  #====================================#
 
  # We select the plot data
  Plot.Data <- Data[Data$Plot==List.Plot[j],]

  # New number of trees for plot
  Levels.IdentTree <- Plot.Data$Tree
  nTree.Plot <- length(Plot.Data$Tree)

  # Computing the distance between trees with dist function
  Dist.Matrix <- as.matrix(dist(cbind(Plot.Data$X,Plot.Data$Y),diag=TRUE,upper=TRUE))

  # Trees with a distance <= to R.Comp
  Include.Mat <- Dist.Matrix
  Include.Mat[Include.Mat<=R.Comp] <- 1
  Include.Mat[Include.Mat>R.Comp] <- 0

  # We do not include the target tree: 0 on diagonal
  for (i in 1:length(Include.Mat[,1])) {
    Include.Mat[i,i] <- 0
  }
  
  # Basal Area vectors: trees with -9999 as diameter have a BA for competition index set to 0 
  for (d in 1:nCens) {
     D.vect <- eval(parse(text=paste("Plot.Data[,5+",d,"]",sep="")))
     D.vect[D.vect==-9999] <- 0
     eval(parse(text=paste("Plot.Data$BA",d," <- pi*(D.vect/200)^2",sep="")))
  }

  # Surface of the neighborhood for each tree
  Sm2 <- surface_comp_index(X.trees=Plot.Data$X,
                            Y.trees=Plot.Data$Y,
                            X.Plot=X.Plot,
                            Y.Plot=Y.Plot,
                            R.Comp=R.Comp)
  Sha <- Sm2/10000 # in ha

  # Competition index for each tree for each date
  for (d in 1:nCens) {
     eval(parse(text=paste("Plot.Data$CI",d," <- Include.Mat%*%Plot.Data$BA",d,"/Sha",sep="")))
  }

  # We replace the DBH at the date of death (0) by 9999
  for (d in 2:nCens) {
    eval(parse(text=paste("Plot.Data[Plot.Data$CensusDeath==",d,",5+",d,"] <- 9999",
                 sep="")))
  }

  # Data.Long
  Data.Long <- data.frame(matrix(NA,ncol=12,nrow=nCens*nTree.Plot))
  colnames(Data.Long) <- c("Plot","Tree","X","Y","Sp","G_tp1","D_t","D_tp1","C_t","t","status_tp1","interval_tp1")
  Data.Long$Plot <- rep(Plot.Data$Plot,each=nCens)
  Data.Long$Tree <- rep(Plot.Data$Tree,each=nCens)
  Data.Long$X <- rep(Plot.Data$X,each=nCens)
  Data.Long$Y <- rep(Plot.Data$Y,each=nCens)
  Data.Long$Sp <- rep(Plot.Data$Species,each=nCens)
  Data.Long$t <- rep(c(1:nCens),nTree.Plot)

  # Filling in Data.Long
  # D_t
  Data.Long$D_t <- c(t(Plot.Data[,c((5+1):(5+nCens))]))
  # D_t+1
  Data.Long$D_tp1 <- c(rbind(t(Plot.Data[,c((5+2):(5+nCens))]),rep(0,nTree.Plot)))
  # C_t
  Data.Long$C_t <- c(t(Plot.Data[,c((7+2*nCens+1):(7+2*nCens+nCens))]))
  # inteval_tp1
  Data.Long$interval_tp1 <- rep(c(Interval,NA),nTree.Plot)
  # G_tp1
  Data.Long$G_tp1 <- (Data.Long$D_tp1-Data.Long$D_t)*10/Data.Long$interval_tp1 # Growth in mm.yr-1
  Data.Long$G_tp1[Data.Long$D_tp1==9999] <- 9999
  # Tree status (0: alive, 1:dead, 2:recruit, not available at t or tp1 => NA)
  Data.Long$status_tp1[(Data.Long$D_t>0&Data.Long$D_tp1>0&
                    Data.Long$D_t!=9999&Data.Long$D_tp1!=9999)] <- 0
  Data.Long$status_tp1[Data.Long$D_tp1==9999] <- 1
  Data.Long$status_tp1[Data.Long$D_t==0&Data.Long$D_tp1>0&Data.Long$D_tp!=9999] <- 2
  Data.Long$status_tp1[Data.Long$D_t==-9999|Data.Long$D_tp1==-9999] <- NA
  # Removing trees without status
  Data.Clean <- Data.Long[!is.na(Data.Long$status_tp1),]

  # Concatenating plot data in Data.pretreat.g.m
  Data.pretreat.g.m <- rbind(Data.pretreat.g.m,Data.Clean)

  #====================================#
  # Preparing Recruitment data
  #====================================#

  # Number of cells on each plot depending on plot size (X.Plot,Y.Plot)
  nXCell <- floor(X.Plot/L.Cell)
  nYCell <- floor(Y.Plot/L.Cell)
  nCell <- nXCell*nYCell
  Lx.LastCell <- L.Cell+X.Plot%%L.Cell
  Ly.LastCell <- L.Cell+Y.Plot%%L.Cell

  # Number of species
  Levels.Species <- sort(unique(as.factor(as.character(Plot.Data$Species))))
  nSpecies <- length(Levels.Species)
  IdentSpecies <- as.numeric(as.factor(as.character(Plot.Data$Species)))

  # Data.Cell
  Data.Cell <- data.frame(matrix(NA,ncol=7+nCens,nrow=nCell))
  names.BA <- paste("BA.Cell.",c(1:nCens),sep="")
  colnames(Data.Cell) <- c("Plot","Cell","X.Cell.min","X.Cell.max","Y.Cell.min","Y.Cell.max","S.Cell",names.BA)
  Data.Cell$Plot <- j
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
  
  # Trees on each cell
  TreesOnCell <- vector()
  for (i in 1:nTree.Plot) {
    TreesOnCell[i] <- min(which(Plot.Data$X[i]>=Data.Cell$X.Cell.min&
                            Plot.Data$X[i]<=Data.Cell$X.Cell.max&
                            Plot.Data$Y[i]>=Data.Cell$Y.Cell.min&
                            Plot.Data$Y[i]<=Data.Cell$Y.Cell.max)) # min allocates the tree in a frontier to the
                                                                   # lower cell number
  }

  #===== Total basal area (m2.ha-1) on each cell at each date
  D.List <- as.matrix(Plot.Data[,(5+1):(5+nCens)])
  D.List[D.List==9999|D.List==-9999] <- 0
  BAr.List <- pi*(D.List/200)^2
  for (c in 1:nCell) {
    Data.Cell[c,(7+1):(7+nCens)] <- c(as.integer(TreesOnCell==c)%*%BAr.List)/(Data.Cell$S.Cell[c]/10000) # Surface must be in ha
  }
  
  #===== Data.recruit
  Data.r.Plot <- as.data.frame(matrix(NA,nrow=nCell*nSpecies*(nCens-1),ncol=9))
  names(Data.r.Plot) <- c("Plot","Cell","Sp","tp1","C_t","BAsp_t","R_tp1","interval_tp1","SCell")
  Data.r.Plot$Plot <- rep(j,nCell*nSpecies*(nCens-1))
  Data.r.Plot$Cell <- rep(c(1:nCell),each=nSpecies*(nCens-1))
  Data.r.Plot$Sp <- rep(rep(Levels.Species,each=nCens-1),nCell)
  Data.r.Plot$tp1 <- rep(c(2:nCens),nSpecies*nCell)
  Data.r.Plot$interval_tp1 <- rep(c(Interval),nSpecies*nCell)
  Data.r.Plot$SCell <- rep(Data.Cell$S.Cell,each=nSpecies*(nCens-1))
  #===== C_t
  Mat.C_t <- matrix(ncol=nCell)
  for (d in 1:(nSpecies)) {
    Mat.C_t <- rbind(Mat.C_t,t(Data.Cell[,(7+1):(7+nCens-1)]))
  }
  Data.r.Plot$C_t <- c(Mat.C_t[-c(1),])
  #===== BAsp_t
  Data.r.Plot$BAsp_t <- 0
  #== table of (cells x species)
  CellSpecies <- as.data.frame(matrix(0,nrow=nCell,ncol=nSpecies))
  names(CellSpecies) <- Levels.Species
  for (c in 1:nCell) {
    for (k in 1:nSpecies) {
      if (Levels.Species[k]%in%Plot.Data$Species[TreesOnCell==c]) {
        CellSpecies[c,k] <- 1
      }
    }
  }
  #== BAsp_t when CellSpecies !=0
  for (c in 1:nCell) {
    for (k in 1:nSpecies) {
      if (CellSpecies[c,k]==1) {
        #== Vector of 1/0 indicating the conspecific trees on cell
        Consp.Tree <- as.integer(TreesOnCell==c & Plot.Data$Species==as.character(Levels.Species[k]))
        #== Vector of TRUE/FALSE for the lines of Data.r.Plot
        Select.Lines <- c(Data.r.Plot$Cell==c & Data.r.Plot$Sp==Levels.Species[k])
        Data.r.Plot$BAsp_t[Select.Lines] <- c(Consp.Tree%*%BAr.List[,-c(nCens)])/(Data.Cell$S.Cell[c]/10000)
      }
    }
  }
  #===== R_tp1
  Data.r.Plot$R_tp1 <- 0
  Data.Clean.r <- Data.Clean[Data.Clean$status_tp1==2,]
  Data.Clean.r$Cell <- 0
  for (i in 1:length(Data.Clean.r$Tree)) {
    Data.Clean.r$Cell[i] <- min(which(Data.Clean.r$X[i]>=Data.Cell$X.Cell.min&
                                      Data.Clean.r$X[i]<=Data.Cell$X.Cell.max&
                                      Data.Clean.r$Y[i]>=Data.Cell$Y.Cell.min&
                                      Data.Clean.r$Y[i]<=Data.Cell$Y.Cell.max))
    Data.r.Plot$R_tp1[Data.r.Plot$Cell==Data.Clean.r$Cell[i]&
                      Data.r.Plot$tp1==(Data.Clean.r$t[i]+1)&
                      Data.r.Plot$Sp==as.character(Data.Clean.r$Sp[i])] <-
                        Data.r.Plot$R_tp1[Data.r.Plot$Cell==Data.Clean.r$Cell[i]&
                                          Data.r.Plot$tp1==(Data.Clean.r$t[i]+1)&
                                          Data.r.Plot$Sp==as.character(Data.Clean.r$Sp[i])]+1
  }

  #===== Merging censuses and replacing C_t and BAsp_t by their initial values
  Data.r.Plot.merged <- Data.r.Plot[order(paste(Data.r.Plot$Plot,Data.r.Plot$Cell,Data.r.Plot$Sp)),]
  for (i in 1:length(DateR)) {
    Data.r.Plot.merged$R_tp1[Data.r.Plot.merged$tp1==DateR[i]] <-
      tapply(Data.r.Plot.merged$R_tp1[Data.r.Plot.merged$tp1 %in% yr.List[[i]]],
             paste(Data.r.Plot.merged$Plot,Data.r.Plot.merged$Cell,Data.r.Plot.merged$Sp)[Data.r.Plot.merged$tp1 %in% yr.List[[i]]],
             sum)
    Data.r.Plot.merged$interval_tp1[Data.r.Plot.merged$tp1==DateR[i]] <- IntervalR[i]
    Data.r.Plot.merged$C_t[Data.r.Plot.merged$tp1==DateR[i]] <- Data.r.Plot.merged$C_t[Data.r.Plot.merged$tp1==yr.List[[i]][1]]
    Data.r.Plot.merged$BAsp_t[Data.r.Plot.merged$tp1==DateR[i]] <- Data.r.Plot.merged$BAsp_t[Data.r.Plot.merged$tp1==yr.List[[i]][1]] 
  }
  
  #===== Selecting censuses and concatenating plot data in Data.pretreat.r
  Data.pretreat.r <- rbind(Data.pretreat.r,Data.r.Plot.merged[Data.r.Plot.merged$tp1%in%DateR,])
}

#====================================#
# Final data-sets
#====================================#

# Growth
Data.growth <- Data.pretreat.g.m[!(Data.pretreat.g.m$status_tp1)%in%c(2,1),
                                 c(1,2,3,4,5,7,9,6)]
write.table(Data.growth,file="data_growth.txt",sep="\t",row.names=FALSE)

# Mortality
Data.mortality <- Data.pretreat.g.m[Data.pretreat.g.m$status_tp1!=2,
                                    c(1,2,3,4,5,7,9,11,12)]
write.table(Data.mortality,file="data_mortality.txt",sep="\t",row.names=FALSE)

# Recruitment
Data.recruitment <- Data.pretreat.r[,c(1,2,3,5,6,7,8,9)]
write.table(Data.recruitment,file="data_recruitment.txt",sep="\t",row.names=FALSE)

cat("\nFunction entry_data has finished, please check files:
1. data_growth.txt
2. data_mortality.txt
3. data_recruitment.txt
in the working directory\n\n")
flush.console()

}

