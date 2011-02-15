#!/usr/bin/env Rscript

# Redirect output:
sink("03-eff-large-pool.log")
pdf("03-eff-large-pool.log.pdf")
# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

nr.cycles 	<- 30
nr.effs		<- 1000
nr.replications	<- 10000

expected.bl<-function(eff, nr.cycles ){
	nr.cycles * ( 1 - 1/(mean(eff) + 1) )
}

mean.bl <- numeric()
exp.bl	<- numeric()


for(i in seq(from=0.5,to=1,by=0.1)){

	efficiencies<-rep(i, nr.cycles)

	x<-PCRcoal(
        	initial.size    =1000,
        	sample.size     =2,
        	nr.cycles       =nr.cycles,
		efficiencies    =efficiencies
	);

	branch.lengths<-numeric()

	for (i in 1:nr.replications){
		branch.lengths[i] <- mean(sample.tree(x)$edge.length)
	}

	mean.bl <- c(mean.bl, mean(branch.lengths) )
	exp.bl <- expected.bl(efficiencies, nr.cycles)
}


plot(exp.bl, mean.bl)

reg<-lm(exp.bl ~ eff.seq)

summary(reg)

plot(reg)

