\name{entry_recruitment}
\alias{entry_recruitment}
\title{
A function to estimate the parameters of the recruitment model for each species in the data-set.
}

\description{ A function to estimate the parameters of the recruitment
  model for each species in the data-set. The function prints summaries
  and plots diagnostic graphics.  }

\usage{entry_recruitment(data_recruitment, burnin=5000, mcmc=5000, thin=5, th.sign=5)}

\arguments{
  \item{data_recruitment}{
    data frame which must have the following column names:
    | Cell | Sp | BAsp_t | C_t | R_tp1 | interval_tp1 | SCell |
  }
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
