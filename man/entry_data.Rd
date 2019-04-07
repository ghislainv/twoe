\name{entry_data}
\alias{entry_data}
\title{
Function to transform the permanent plot data-set into right format data-sets for the growth, mortality and recruitment model.
}
\description{
Function to transform the permanent plot data-set into right format data-sets for the growth, mortality and recruitment model.
}
\usage{
entry_data(Data, XY.Plot, R.Comp, L.Cell)
}
\arguments{
  \item{Data}{
Data frame of the permanent plot data-set.
}
  \item{XY.Plot}{
Data frame with the dimensions (m) of each plot: | Plot | X.Plot |
Y.Plot |
}
  \item{R.Comp}{
radius (m) for the competition area around each tree
}
  \item{L.Cell}{
width (m) of the ground cell for recruitment estimation
}
}
\value{
Return three data-sets: data_growth.txt, data_mortality.txt and data_recruitment.txt
}
\author{
Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}
\keyword{Permanent plot data-set}
