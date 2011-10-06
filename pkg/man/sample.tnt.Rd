\name{sample.tnt}
\alias{sample.tnt}
\alias{sample.tnt,PCRcoal-method}

\title{Sample a coalescent tree according to a PCR experiment, return also subsample sizes and size trajectories}

\description{
	This method returns a coalescent tree sampled according to the PCR experiment specified by the 
	parameters stored in the slots of the \code{\link{PCRcoal}} object given as the first argument.
	All of the slots must contain a defined value.

	This method also returns the subsample sizes and size trajectories corresponding to the sampled tree.
}

\usage{
	\S4method{sample.tnt}{PCRcoal}(object)
}

\arguments{
  \item{object}{A \code{\link{PCRcoal}} object, with all the slots defined.}
}

\value{
	A list containing the following elements:
	\describe{
	\item{phylo}{A \code{\link{ape}} \code{phylo} object (rooted).}
	\item{trajectories}{A matrix containing the sampled size trajectories (rows ~ subsamples, columns ~ cycles).}
	\item{subsamples}{A vector of subsample sizes.}
	}
}

\seealso{The \code{\link{PCRcoal}} class. The \code{\link{sample.tree}} and \code{\link{sample.trs}} methods.}


\examples{
	# Construct a PCRcoal object 
	# and specify the experimental conditions:
	sim <-PCRcoal(
        initial.size    =5,
        sample.size     =10,
        nr.cycles       =30,
        efficiencies    =c(rep(0.1, 30))
    	)

	# Sample a coalescent tree with subsample sizes and
    # size trajectories:
	res <- sample.tnt(sim)

	# Print subsample sizes:
	print(res$subsamples)
	
	# Print size trajectories:
	print(res$trajectories)

	# Print tree info:
	print(res$phylo)

	# Plot the tree:
	plot(res$phylo)
}

\keyword{PCR, simulation}
