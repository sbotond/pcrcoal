#!/usr/bin/env Rscript

# This scripts runs a large number of coalescent PCR simulation
# with a single sampled molecule and reports basic statistics about
# the distribution of the length of the simulated branches.
#
# Output:
#   * 01-bl-test.log - summary of the branch length distribution
#   * 01-bl-test.pdf - histogram of simulated branch lengths

# The probability that the simulated lineage "passes through" a newly synthesised
# molecule (which increases branch length by one) in a given cycle is proportional 
# to the number of such molecules present after the cycle.
#
# This can be calculated as eff/(eff + 1), where "eff" is the per-cycle efficiency.
# Hence, the expected branch length under a fixed efficiency across cycles
# is [eff/(eff +1 )] * nr_cycles.

# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

# Redirect output:
sink("01-bl-test.log")
pdf("01-bl-test.pdf")

nr_replications     <- 10000            # Number of simulated PCR experiments
effs                <- c(rep(0.5, 30))  # Per-cycle PCR efficiencies
nr.cycles           <- 30               # Number of PCR cycles

x<-PCRcoal(
        initial.size    =10,    # Initial number of molecules
        sample.size     =1,     # Number of sampled molecules
        nr.cycles       =nr.cycles,    # Number of PCR cycles
	    efficiencies    =effs   # Per-cycle PCR efficiencies
);

branch.lengths<-numeric()

for (i in 1:nr_replications){
	branch.lengths[i] <- sample.tree(x)$edge.length[1]
}

# Calculate expected value of the simulated branch length:
calc_exp_bl<-function(eff, nr_cycles){
   return(nr_cycles * ( eff/(eff +1 ))) 
}

exp_bl <- calc_exp_bl(mean(effs), nr.cycles)

cat("Expected mean branch length:", exp_bl, "\n")
print(summary(branch.lengths))


hist(branch.lengths,nclass=30)
abline(v=exp_bl,col="red")

