#!/usr/bin/env Rscript

# Script for profiling the sample.tree method.
#
# Output: pcrcoal_prof.out, to be analysed with R CMD Rprof.

library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");
source("../PCRcoalSource.R");

x<-PCRcoal(
    initial.size    =1000,
    sample.size     =500,
    nr.cycles       =33,
    eff    =rep(0.5, 33)
);

Rprof("pcrcoal_prof.out")
t<-sample.tree(x);
Rprof(NULL)

system("R CMD Rprof pcrcoal_prof.out")

