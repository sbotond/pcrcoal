#!/usr/bin/env Rscript

# This script runs a large number of coalescent PCR simulations with one sampled molecule.
# The efficiency vectors are sampled from the range [0.5, 1.0].
# It runs several simulation under a sampled efficiency vector and saves the mean simulated
# branch length.
#
# The expected value of the simulated branch lengths is approximated as nr_cycles * [ mean(eff) / (1 + mean(eff))],
# where "eff" is the vector of per-cycle PCR efficiencies (see also 01-bl-test.R).
#
# Output:
#   * 02-bl-eff-test.pdf - scatter plot of expected versus mean realized branch lengths under a given efficiency vector
#   * 02-bl-eff-test.log - linear regression of mean realized vs expected branch lengths

# Redirect output:
sink("02-bl-eff-test.log")
pdf("02-bl-eff-test.pdf")

# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

# Numberf of PCR cycles:
nr.cycles 	<- 30
# Number of sampled efficiency vectors:
nr.effs		<- 1000
# Number of "trees" sampled for a given efficiency vector:
nr.replications	<- 5000

# Function to calculate the expected branch length:
expected.bl<-function(eff, nr.cycles ){
	nr.cycles * ( eff/(eff + 1) )
}

# Function to sample a vector of efficiencies having a mean between 0.5 and 1.0:
sample.effs <- function(nr.cycles){
	ru <- runif(1,min=0.5,max=1)
        effs <- rbeta(nr.cycles, shape1=1,shape2=((1-ru)/ru) )
	effs <- (effs * 0.5) + 0.5
        return(effs)
}

mean.bl <- numeric()
exp.bl	<- numeric()

for(i in 1:nr.effs){

	# Sample efficiencies:
	efficiencies<-sample.effs(nr.cycles)

	# Construct the simulation object:
	x<-PCRcoal(
        	initial.size    =1,
        	sample.size     =1,
        	nr.cycles       =nr.cycles,
		efficiencies    =efficiencies
	);

	branch.lengths<-numeric()
	
	# Sample nr.replications branch lengths:
	for (i in 1:nr.replications){
		branch.lengths[i] <- sample.tree(x)$edge.length[1]
	}

	# Update the expected and observed branch length vectors:
	mean.bl <- c( mean.bl, mean(branch.lengths) )
	exp.bl <- c(exp.bl, expected.bl(mean(efficiencies), nr.cycles))
}

# Plot expected vs. observed branch lenghts:
plot(mean.bl, exp.bl)

# Regression through (0,0):
reg<-lm(mean.bl ~ -1 + exp.bl)

summary(reg)

plot(reg)


