\name{surface_comp_index}
\alias{surface_comp_index}
\title{Function to compute the competition area (m^2) around each target tree.
}
\description{
Function to compute the competition area (m^2) around each target tree. As some trees are close to the plot edges, the competition area is changing from tree to tree. The competition area is equal to the intersection between a circle having the target tree as center and the rectangle delimited by the plot edges.
}
\usage{
surface_comp_index(X.trees, Y.trees, X.Plot, Y.Plot, R.Comp)
}
\arguments{
  \item{X.trees}{
vector of X coordinates (m) of each tree
}
  \item{Y.trees}{
vector of Y coordinates (m) of each tree
}
  \item{X.Plot}{
size of the plot (m) on the X axis
}
  \item{Y.Plot}{
size of the plot (m) on the Y axis
}
  \item{R.Comp}{
radius (m) for the competition area around each tree
}
}
\value{
Return a vector of competition area (m^2) for each tree
}
\author{
Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}
\examples{
\dontrun{
X.trees <- runif(50,0,100) # X coordinates of the trees
Y.trees <- runif(50,0,100)  # Y coordinates of the trees
X.Plot <- 100 # Plot is a 1ha plot with edge=100m
Y.Plot <- 100
R.Comp <-10 # Competition on a 10m radius around each tree
S_vect <- surface_comp_index(X.trees, Y.trees, X.Plot, Y.Plot, R.Comp)
}
}
\keyword{competition index}
