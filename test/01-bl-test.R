#!/usr/bin/env Rscript

# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

# Redirect output:
sink("01-bl-test.log")
pdf("01-bl-test.pdf")

x<-PCRcoal(
        initial.size    =1,
        sample.size     =1,
        nr.cycles       =30,
	efficiencies    =c(rep(0.5, 30))
);

branch.lengths<-numeric()

for (i in 1:10000){
	branch.lengths[i] <- sample.tree(x)$edge.length[1]
}

print(summary(branch.lengths))

hist(branch.lengths,nclass=30)


