\name{sample.tree}
\alias{sample.tree}
\alias{sample.tree,PCRcoal-method}

\title{Sample a coalescent tree according to a PCR experiment}

\description{
	This method returns a coalescent tree sampled according to the PCR experiment specified by the 
	parameters stored in the slots of the \code{\link{PCRcoal}} object given as the first argument.
	All of the slots must contain a defined value.

	Use this method if you are interested only in the sampled tree.
}

\usage{
	\S4method{sample.tree}{PCRcoal}(object)
}

\arguments{
  \item{object}{A \code{\link{PCRcoal}} object, with all the slots defined.}
}

\value{
	An \code{\link{ape}} \code{phylo} object (rooted).
}

\seealso{The \code{\link{PCRcoal}} class. The \code{\link{sample.tnt}} and \code{\link{sample.trs}} methods.}

\examples{
	# Construct a PCRcoal object 
	# and specify the experimental conditions:
	sim <-PCRcoal(
        initial.size    =5,
        sample.size     =10,
        nr.cycles       =30,
        efficiencies    =c(rep(0.1, 30))
    	)

	# Sample a coalescent tree:
	tree <- sample.tree(sim)

	# Print tree info:
	print(tree)

	# Plot the tree:
	plot(tree)
	
}

\keyword{PCR, simulation}
