\name{PCRcoal-class}
\docType{class}
\alias{PCRcoal}
\alias{PCRcoal-class}
\alias{PCRcoal,missing-method}


\title{The PCRcoal class}

\description{

PCRcoal is a class implementing the coalescent approach to PCR simulation developed by Weiss and von Haeseler (NAR, 1997).

}

\note{
\itemize{
	\item The original method of Weiss and von Haeseler produces a tree containing a multifurcating root node connected to the subtrees
	corresponding to the different subsamples. In order to be able to store the output in \code{phylo} objects, \code{PCRcoal} randomly breaks the polytomy by "coalescing" the subsample trees with zero branch lengths. Some branch length might be non-zero, in the case a subsample coalesced into one molecule which accumulated branch length ("replication count") in the previous cycles.

	\item If the cycle number and/or the efficiencies are small, it can happen that the number of the molecules after the amplification is less than the specified sample size. \code{PCRcoal} will try sampling size trajectories \code{max.tries} times (100 by default) to get enough amplified molecules and then aborts the simulation.
}
}

\references{
	Weiss, G, von Haeseler, A (1997) A coalescent approach to the polymerase chain reaction. NAR 25(15): 3082-3087 \url{http://dx.doi.org/10.1093/nar/25.15.3082} 
}

\section{Usage}{
	\code{
	PCRcoal(
        initial.size,
		nr.cycles,
		efficiencies,
		sample.size,
		max.tries=100
	)
	}

	Returns a \code{PCRcoal} S4 object.
}

\section{Slots}{
  \describe{
    	\item{\code{initial.size}:}{Number of template molecules.}
    	\item{\code{nr.cycles}:}{Number of PCR cycles.}
    	\item{\code{efficiencies}:}{A vector of per-cycle PCR efficiencies. It is recycled if its length is smaller than \code{nr.cycles}.}
    	\item{\code{sample.size}:}{Number of molecules sampled after amplification.}
    	\item{\code{max.tries}:}{Maximum number of tries when sampling a size trajectory with a large enough final size.}
    }
}

\section{Methods}{
  \describe{
    	\item{\code{sample.tree}:}{Sample a tree, return as an \code{\link{ape}} \code{phylo} object.}
    	\item{\code{sample.tnt}:}{Return the sampled tree, subsample sizes and size trajectories in a list.}
    	\item{\code{sample.trs}:}{Return the sampled subsample sizes and size trajectories in a list.}
    }
}

\seealso{
  The \code{\link{ape}} and \code{\link{phylosim}} packages.}
\author{Botond Sipos}
\examples{
	# Construct a PCRcoal object 
        # and specify the experimental conditions:
        sim <-PCRcoal(
        initial.size    =5,
        sample.size     =10,
        nr.cycles       =20,
        efficiencies    =c(rep(0.1, 30))
        )

	##
	## Method: sample.tree
	##

	 # Sample a coalescent tree:
        tree <- sample.tree(sim)

        # Print tree info:
        print(tree)

        # Plot the tree:
        plot(tree)
	
	##
	## Method: sample.tnt
	##

	# Sample a coalescent tree with subsample sizes and
	# size trajectories:
        res <- sample.tnt(sim)

        # Print subsample sizes:
        print(res$subsamples)

        # Print size trajectories:
        print(res$trajectories)

        # Print tree info:
        print(res$phylo)
	
	##
	## Method: sample.trs
	##

	# Sample subsample sizes and size trajectories:
        res <- sample.trs(sim)

        # Print subsample sizes:
        print(res$subsamples)

        # Print size trajectories:
        print(res$trajectories)

	##
	## Simulating mutations using PhyloSim
	##

	# Sample tree:
	tree <- sample.tree(sim)

	# Load the phylosim package
	library(phylosim)
	# Construct a PhyloSim object, set up the root
	# sequence and the substitution process:	
	psim <- PhyloSim(
		phylo = tree,
		root.seq=sampleStates( NucleotideSequence(length=10,processes=list(list(JC69())) ) )
	)
	# Scale the tree:	
	scaleTree(psim,0.1)

	# Simulate mutations:
	Simulate(psim)

	# Plot the resulting alignment:
	plot(psim,plot.ancestors=FALSE)
}

\keyword{classes, coalescent, PCR, simulation}
