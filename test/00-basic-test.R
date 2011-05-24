#!/usr/bin/env Rscript

# This script will simulate coalescent PCR trees under a 
# range of parameters. It takes a long time to finish.
#
# Output: 00-basic-test.log

# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Redirect output:
sink("00-basic-test.log")

# Load package source:
source("../PCRcoalSource.R");

init.sizes      <-1:500 # Initial number of molecules.
sample.sizes    <-2:500 # Number of sampled molecules after PCR.
nr.cycles       <-10:25 # Number of cycles.
effs            <-seq(from=0.2,to=1.0,by=0.1) # Per-cycle efficiencies.

for(i.sz    in init.sizes){
for(s.sz    in sample.sizes){
for(nr.c    in nr.cycles){
for(ef      in effs){

    if(s.sz > nr.c){ next }

    cat("Initial size:\t", i.sz,"\n")
    cat("Sample size:\t", s.sz,"\n")
    cat("Nr of cycles:\t", nr.c,"\n")
    cat("Efficiency:\t", ef,"\n")
    cat("\n")

    x<-PCRcoal(
        initial.size    =i.sz,
        sample.size     =s.sz,
        nr.cycles       =nr.c,
        efficiencies    =c(rep(ef, nr.c))
    );

    t<-sample.tree(x)

}
}
}
}


