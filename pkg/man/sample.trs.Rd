\name{sample.trs}
\alias{sample.trs}
\alias{sample.trs,PCRcoal-method}

\title{Sample subsample sizes and size trajectories according to a PCR experiment}

\description{
	This method also returns the subsample sizes and size trajectories sampled according to the PCR experiment specified by the
        parameters stored in the slots of the \code{\link{PCRcoal}} object given as the first argument.

	This method does not sample the tree, which makes it faster than the \code{\link{sample.tree}} method.
	Use this method if you are interested only in sampling subsamples and size trajectories.
}

\usage{
	\S4method{sample.trs}{PCRcoal}(object)
}

\arguments{
  \item{object}{A \code{\link{PCRcoal}} object, with all the slots defined.}
}

\value{
	A list containing the following elements:
	\describe{
	\item{trajectories}{A matrix containing the sampled size trajectories (rows ~ subsamples, columns ~ cycles).}
	\item{subsamples}{A vector of subsample sizes.}
	}
}

\seealso{The \code{\link{PCRcoal}} class. The \code{\link{sample.tree}} and \code{\link{sample.tnt}} methods.}


\examples{
	# Construct a PCRcoal object 
	# and specify the experimental conditions:
	sim <-PCRcoal(
        initial.size    =5,
        sample.size     =10,
        nr.cycles       =30,
        efficiencies    =c(rep(0.1, 30))
    	)

	# Sample subsample sizes and size trajectories:
	res <- sample.trs(sim)

	# Print subsample sizes:
	print(res$subsamples)

	# Print size trajectories:
	print(res$trajectories)

}

\keyword{PCR, simulation}
