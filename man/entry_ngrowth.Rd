\name{entry_ngrowth}
\alias{entry_ngrowth}
\title{
  A function to estimate the parameters of the growth model for each species in the data-set.
}

\description{ A function to estimate the parameters of the growth model
  for each species in the data-set. The function prints summaries and
  plots diagnostic graphics.}

\usage{entry_ngrowth(data_growth, a.sd=0.927, b.sd=0.0038, burnin=1000,
mcmc=1000, thin=1, th.sign=5)}

\arguments{
  \item{data_growth}{
    data frame which must have the following column names: | Tree | Sp | D_t | C_t | G_tp1 |
  }

  \item{a.sd}{A positive scalar for the intercept of the function
  describing the standard deviation for the observation error
  process. Default to 0.927 (see Ruger et al. 2005).}

  \item{b.sd}{A positive scalar for the slope of the function describing
  the standard deviation for the observation error process. Default to
  0.0038 (see Ruger et al. 2005).}
  
  \item{burnin}{The number of burnin iterations for the sampler.}
  
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
  
  \item{thin}{The thinning interval used in the simulation. The number of
    mcmc iterations must be divisible by this value.}

  \item{th.sign}{Threshold (in \%) used to test parameter significance:
    is zero included in the interval defined by the threshold? Default
    to 5.}
}
\value{
The function prints summaries and plots diagnostic graphics.
}
\author{
Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}
\keyword{Gibbs sampler}
