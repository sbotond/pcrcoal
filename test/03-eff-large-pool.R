#!/usr/bin/env Rscript

# Redirect output:
sink("03-eff-large-pool.log")
pdf("03-eff-large-pool.pdf")

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
	nr.cycles * ( 1 - 1/(eff + 1) )
}

# Function to sample a vector of efficiencies having a mean between 0.5 and 1.0:
sample.effs <- function(nr.cycles){
	ru <- runif(1,min=0,max=1)
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
        	initial.size    =3000,
        	sample.size     =3,
        	nr.cycles       =nr.cycles,
		efficiencies    =efficiencies
	);

	branch.lengths<-numeric()
	
	# Sample nr.replications branch lengths, save the mean of the non-zero branch lengths:
	for (i in 1:nr.replications){
		edges<-sample.tree(x)$edge.length
		branch.lengths[i] <- mean(edges[ which(edges > 0) ])
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

