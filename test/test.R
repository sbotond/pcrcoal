#!/usr/bin/env Rscript

# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

init.sizes      <-1:500
sample.sizes    <-2:500
nr.cycles       <-10:25
effs            <-seq(from=0.2,to=1.0,by=0.1)

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


