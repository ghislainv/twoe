\name{exit_simu}
\alias{exit_simu}
\title{
Function to simulate tropical forest dynamics.
}
\description{
Function to simulate tropical forest dynamics. It uses the species
parameters obtain by the "exit_..." functions. Function also plots the
evolution of the basal area of each selected plot for simulations
(exit_BA.pdf) and the evolution of the basal area of each species on
each plot (exit_biodiversity.pdf).
}
\usage{
exit_simu(Data, XY.Plot, D.Recruitment, R.Comp, L.Cell, Plot.Sim,
Year.Sim, Step.Sim)
}
\arguments{
  \item{Data}{
Data frame of the permanent plot data-set.
}
  \item{XY.Plot}{
Data frame with the dimensions (m) of each plot: |Plot|X.Plot|Y.Plot|
}
  \item{D.Recruitment}{
recruitment diameter (cm)
}
  \item{R.Comp}{
radius (m) for the competition area around each tree
}
  \item{L.Cell}{
width (m) of the ground cell for recruitment estimation
}
  \item{Plot.Sim}{
vector of plots name choosen for simulations
}
  \item{Year.Sim}{
number of years of simulation from the first census 
}
  \item{Step.Sim}{
time step for which the simulation are kept in memory  
}
}
\value{
Return a data-set with the monitoring of each tree from first census and
plot the evolution of the basal area of the plot with time.
}
\author{
Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}
\keyword{forest dynamics simulation}
